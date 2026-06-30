# =============================================================================
# modeling_functions.R
#
# Shared modeling helpers used by NCSS_models.qmd and Perry_predict.qmd.
# Source this file from both scripts.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(tibble); library(forcats)
  library(tidymodels); library(recipes); library(yardstick)
  library(uwot); library(dbscan); library(rlang); library(workflows)
  library(tune); library(dials); library(parsnip); library(xgboost)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# -----------------------------------------------------------------------------
# Columns we drop from predictors (every taxonomy label is leakage)
# -----------------------------------------------------------------------------
TAXONOMY_DROP_COLS <- c(
  "taxonname_clean", "taxonname_key", "taxclname",
  "taxorder", "taxsuborder", "taxgrtgroup",
  "taxsubgrp", "taxsubgrp_mod1", "taxsubgrp_mod2", "taxsubgrp_mod3",
  "taxsubgrp_mod_group", "taxgrtgroup_mod_group", "taxsuborder_mod_group",
  "taxpartsize"
)

# -----------------------------------------------------------------------------
# drop_bad_predictors: tidies up a train/test pair before fitting
# -----------------------------------------------------------------------------
drop_bad_predictors <- function(train_df, test_df, outcome,
                                id_col = "peiid", drop_cols = NULL) {
  train_df <- train_df %>% select(-any_of(drop_cols))
  test_df  <- test_df  %>% select(-any_of(drop_cols))
  
  train_df <- train_df %>% mutate(across(where(is.logical), as.integer))
  test_df  <- test_df  %>% mutate(across(where(is.logical), as.integer))
  
  char_preds <- setdiff(
    names(train_df)[vapply(train_df, is.character, logical(1))],
    c(id_col, outcome)
  )
  if (length(char_preds) > 0) {
    train_df <- train_df %>% select(-all_of(char_preds))
    test_df  <- test_df  %>% select(-all_of(char_preds))
  }
  
  pred_names <- setdiff(names(train_df), c(id_col, outcome))
  one_val <- pred_names[
    vapply(train_df[pred_names],
           function(x) length(unique(x[!is.na(x)])) <= 1, logical(1))
  ]
  if (length(one_val) > 0) {
    train_df <- train_df %>% select(-all_of(one_val))
    test_df  <- test_df  %>% select(-all_of(one_val))
  }
  
  list(train = train_df, test = test_df)
}

# -----------------------------------------------------------------------------
# prep_outcome: filter rare classes, factorize outcome, drop singletons
# -----------------------------------------------------------------------------
prep_outcome <- function(df, outcome, min_n = 10) {
  df %>%
    filter(!is.na(.data[[outcome]])) %>%
    mutate(!!outcome := factor(.data[[outcome]])) %>%
    add_count(.data[[outcome]], name = ".class_n") %>%
    filter(.class_n >= min_n) %>%
    select(-.class_n) %>%
    droplevels()
}

# -----------------------------------------------------------------------------
# Recipe builders
# -----------------------------------------------------------------------------
make_recipe_for_engine <- function(train_df, outcome) {
  recipe(stats::as.formula(paste(outcome, "~ .")), data = train_df) %>%
    update_role(peiid, new_role = "id") %>%
    step_zv(all_predictors()) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_normalize(all_numeric_predictors())
}

# -----------------------------------------------------------------------------
# bake_xy: prep recipe and produce numeric matrices for UMAP
# -----------------------------------------------------------------------------
bake_xy <- function(train_df, test_df, outcome) {
  rec <- make_recipe_for_engine(train_df, outcome)
  prep_rec <- prep(rec, training = train_df)
  
  bake_one <- function(df) {
    bake(prep_rec, new_data = df) %>%
      select(-any_of(c(outcome, "peiid"))) %>%
      mutate(across(where(is.logical), as.double)) %>%
      select(where(is.numeric)) %>%
      mutate(across(everything(), ~ ifelse(is.finite(.x), .x, 0)))
  }
  
  list(
    prep_rec = prep_rec,
    x_train  = bake_one(train_df),
    x_test   = bake_one(test_df),
    y_train  = train_df[[outcome]],
    y_test   = test_df[[outcome]]
  )
}

# -----------------------------------------------------------------------------
# UMAP + HDBSCAN feature builder (training / scoring symmetric)
# -----------------------------------------------------------------------------
fit_umap_hdb <- function(x_train, x_test, y_train, y_test, outcome,
                         n_components = 10, n_neighbors = 15,
                         min_dist = 0.1, metric = "euclidean",
                         minPts = 10, seed = 2) {
  
  set.seed(seed)
  umap_fit <- uwot::umap(
    as.matrix(x_train),
    n_components = n_components,
    n_neighbors  = n_neighbors,
    min_dist     = min_dist,
    metric       = metric,
    ret_model    = TRUE
  )
  
  z_train <- as.data.frame(umap_fit$embedding)
  names(z_train) <- paste0("U", seq_len(ncol(z_train)))
  z_test <- as.data.frame(uwot::umap_transform(as.matrix(x_test), umap_fit))
  names(z_test) <- names(z_train)
  
  hdb_fit <- dbscan::hdbscan(as.matrix(z_train), minPts = minPts)
  
  train_mat <- as.matrix(z_train)
  test_mat  <- as.matrix(z_test)
  test_cluster <- apply(test_mat, 1, function(z) {
    d <- rowSums((train_mat - matrix(z, nrow(train_mat), ncol(train_mat),
                                     byrow = TRUE))^2)
    hdb_fit$cluster[which.min(d)]
  })
  levs <- sort(unique(c(hdb_fit$cluster, test_cluster)))
  
  train_dat <- bind_cols(
    z_train,
    setNames(list(y_train), outcome),
    tibble(
      hdb_cluster = factor(hdb_fit$cluster, levels = levs),
      hdb_noise   = factor(if_else(hdb_fit$cluster == 0, "noise", "cluster"),
                           levels = c("cluster", "noise"))
    )
  )
  test_dat <- bind_cols(
    z_test,
    setNames(list(y_test), outcome),
    tibble(
      hdb_cluster = factor(test_cluster, levels = levs),
      hdb_noise   = factor(if_else(test_cluster == 0, "noise", "cluster"),
                           levels = c("cluster", "noise"))
    )
  )
  
  list(train = train_dat, test = test_dat,
       umap_fit = umap_fit, hdb_fit = hdb_fit, z_train = z_train)
}

# -----------------------------------------------------------------------------
# Model specs
# -----------------------------------------------------------------------------
make_rf_spec <- function() {
  rand_forest(trees = 1000, mtry = tune(), min_n = tune()) %>% #RF performance plateaus quickly with trees; 1000 is well past the plateau for any reasonable problem size.
    set_engine("ranger", probability = TRUE, importance = "permutation") %>%
    set_mode("classification")
}

make_xgb_spec <- function() {
  boost_tree(
    trees          = tune(), tree_depth     = tune(),
    learn_rate     = tune(), loss_reduction = tune(),
    min_n          = tune(), sample_size    = tune(),
    mtry           = tune(), stop_iter      = 25
  ) %>%
    set_engine("xgboost") %>%
    set_mode("classification")
}

make_grid <- function(engine, n_predictors, size = 12) {
  if (engine == "rf") {
    grid_regular(
      mtry(range  = c(2L, min(20L, n_predictors))),
      min_n(range = c(2L, 25L)),
      levels = ceiling(sqrt(size))
    )
  } else {
    params <- parameters(
      trees(range          = c(200L, 1500L)),
      tree_depth(range     = c(2L, 8L)),
      learn_rate(range     = c(-3, -0.5)),
      loss_reduction(range = c(-2, 0)),
      min_n(range          = c(2L, 15L)),
      sample_prop(range    = c(0.6, 1.0)),
      finalize(mtry(range  = c(2L, min(30L, n_predictors))),
               data.frame(matrix(0, nrow = 1, ncol = n_predictors)))
    )
    grid_latin_hypercube(params, size = size)
  }
}

# -----------------------------------------------------------------------------
# Inner CV tuning + final fit on a single train/test pair
# -----------------------------------------------------------------------------
fit_one_model <- function(train_df, test_df, outcome,
                          engine = c("rf", "xgb"),
                          space  = c("plain", "umap_hdb"),
                          inner_v = 5, grid_size = 12, seed = 1, #grid size = 12 due to computational budget
                          umap_settings = list()) {
  
  engine <- match.arg(engine)
  space  <- match.arg(space)
  set.seed(seed)
  
  # Standardize outcome name for downstream consistency
  train_df <- train_df %>% rename(.y = all_of(outcome))
  test_df  <- test_df  %>% rename(.y = all_of(outcome))
  
  if (space == "plain") {
    train_used <- train_df
    test_used  <- test_df
    rec <- recipe(.y ~ ., data = train_used) %>%
      update_role(peiid, new_role = "id") %>%
      step_mutate_at(where(is.logical), fn = as.numeric) %>%   # <-- ADD THIS LINE
      step_zv(all_predictors()) %>%
      step_impute_median(all_numeric_predictors()) %>%
      step_normalize(all_numeric_predictors())
    prep_rec <- prep(rec, training = train_used)
    n_pred <- ncol(bake(prep_rec, new_data = train_used)) - 2  # -outcome, -id
    
    spec <- if (engine == "rf") make_rf_spec() else make_xgb_spec()
    wf <- workflow() %>% add_model(spec) %>% add_recipe(rec)
    
    grid <- make_grid(engine, n_pred, size = grid_size)
    
    v <- min(inner_v, min(count(train_used, .y)$n))
    folds <- vfold_cv(train_used, v = v, strata = .y)
    
    tuned <- tune_grid(
      wf, resamples = folds, grid = grid,
      metrics = metric_set(kap, bal_accuracy, accuracy)
    )
    best <- select_best(tuned, metric = "kap")
    final_fit <- finalize_workflow(wf, best) %>% fit(train_used)
    
    pred <- predict(final_fit, test_used, type = "prob") %>%
      bind_cols(predict(final_fit, test_used, type = "class")) %>%
      bind_cols(test_used %>% select(.y)) %>%
      rename(!!outcome := .y)
    
    return(list(
      engine = engine, space = space, outcome = outcome,
      tuned = tuned, best = best, final_fit = final_fit,
      prep_rec = prep_rec, pred = pred,
      train_raw = train_used %>% rename(!!outcome := .y),
      umap_fit = NULL, hdb_fit = NULL, z_train = NULL
    ))
  }
  
  # ---- UMAP+HDBSCAN space ----
  uset <- modifyList(
    list(n_components = 10, n_neighbors = 15,
         min_dist = 0.1, metric = "euclidean", minPts = 10),
    umap_settings
  )
  
  # bake into numeric features for UMAP
  rec_pre <- recipe(.y ~ ., data = train_df) %>%
    update_role(peiid, new_role = "id") %>%
    step_mutate_at(where(is.logical), fn = as.numeric) %>%   # <-- ADD THIS LINE
    step_zv(all_predictors()) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_normalize(all_numeric_predictors())
  prep_rec <- prep(rec_pre, training = train_df)
  
  bake_num <- function(df) {
    bake(prep_rec, new_data = df) %>%
      select(-any_of(c(".y", "peiid"))) %>%
      mutate(across(where(is.logical), as.double)) %>%
      select(where(is.numeric)) %>%
      mutate(across(everything(), ~ ifelse(is.finite(.x), .x, 0)))
  }
  x_train <- bake_num(train_df)
  x_test  <- bake_num(test_df)
  
  uh <- fit_umap_hdb(
    x_train, x_test, train_df$.y, test_df$.y, outcome = ".y",
    n_components = uset$n_components, n_neighbors = uset$n_neighbors,
    min_dist = uset$min_dist, metric = uset$metric,
    minPts = uset$minPts, seed = seed + 1
  )
  
  spec <- if (engine == "rf") make_rf_spec() else make_xgb_spec()
  wf <- workflow() %>% add_model(spec) %>% add_formula(.y ~ .)
  
  v <- min(inner_v, min(count(uh$train, .y)$n))
  folds <- vfold_cv(uh$train, v = v, strata = .y)
  
  n_pred <- ncol(uh$train) - 1
  grid <- make_grid(engine, n_pred, size = grid_size)
  
  tuned <- tune_grid(
    wf, resamples = folds, grid = grid,
    metrics = metric_set(kap, bal_accuracy, accuracy)
  )
  best <- select_best(tuned, metric = "kap")
  final_fit <- finalize_workflow(wf, best) %>% fit(uh$train)
  
  pred <- predict(final_fit, uh$test, type = "prob") %>%
    bind_cols(predict(final_fit, uh$test, type = "class")) %>%
    bind_cols(uh$test %>% select(.y)) %>%
    rename(!!outcome := .y)
  
  list(
    engine = engine, space = space, outcome = outcome,
    tuned = tuned, best = best, final_fit = final_fit,
    prep_rec = prep_rec, pred = pred,
    train_raw = train_df %>% rename(!!outcome := .y),
    umap_fit = uh$umap_fit, hdb_fit = uh$hdb_fit, z_train = uh$z_train
  )
}

# -----------------------------------------------------------------------------
# Outer-loop nested CV: returns one fold-level metrics row per resample
# -----------------------------------------------------------------------------
nested_cv_one_variant <- function(df, outcome,
                                  engine, space,
                                  outer_v = 5, inner_v = 5,
                                  grid_size = 12, seed = 1,
                                  umap_settings = list()) {
  
  set.seed(seed)
  outer_folds <- vfold_cv(df, v = outer_v, strata = !!sym(outcome))
  
  fold_results <- map_dfr(seq_len(nrow(outer_folds)), function(i) {
    train_i <- analysis(outer_folds$splits[[i]])
    test_i  <- assessment(outer_folds$splits[[i]])
    
    # Re-droplevels to handle any class lost in the train fold
    train_i <- train_i %>% droplevels()
    test_i  <- test_i  %>% mutate(!!outcome := factor(
      .data[[outcome]], levels = levels(train_i[[outcome]])
    )) %>% filter(!is.na(.data[[outcome]]))
    
    if (nrow(test_i) == 0 || nlevels(droplevels(train_i[[outcome]])) < 2) {
      return(tibble(fold = i, accuracy = NA_real_, kap = NA_real_,
                    bal_accuracy = NA_real_, f_meas_macro = NA_real_,
                    n_test = nrow(test_i)))
    }
    
    res <- tryCatch(
      fit_one_model(
        train_i, test_i, outcome = outcome,
        engine = engine, space = space,
        inner_v = inner_v, grid_size = grid_size,
        seed = seed + i, umap_settings = umap_settings
      ),
      error = function(e) NULL
    )
    if (is.null(res)) {
      return(tibble(fold = i, accuracy = NA_real_, kap = NA_real_,
                    bal_accuracy = NA_real_, f_meas_macro = NA_real_,
                    n_test = nrow(test_i)))
    }
    
    pr <- res$pred %>% rename(.truth = all_of(outcome))
    tibble(
      fold = i,
      n_test = nrow(pr),
      accuracy     = accuracy(pr,     truth = .truth, estimate = .pred_class)$.estimate,
      kap          = kap(pr,          truth = .truth, estimate = .pred_class)$.estimate,
      bal_accuracy = bal_accuracy(pr, truth = .truth, estimate = .pred_class)$.estimate,
      f_meas_macro = f_meas(pr,       truth = .truth, estimate = .pred_class,
                            estimator = "macro")$.estimate
    )
  })
  
  fold_results %>%
    mutate(engine = engine, space = space, outcome = outcome)
}

# -----------------------------------------------------------------------------
# Bundle constructors (Perry-deployment compatible)
# -----------------------------------------------------------------------------
build_plain_bundle <- function(name, outcome, train_raw, prep_rec, final_fit) {
  x_cols <- bake(prep_rec, new_data = train_raw) %>%
    select(-any_of(c(outcome, "peiid"))) %>%
    mutate(across(where(is.logical), as.double)) %>%
    select(where(is.numeric)) %>%
    names()
  
  list(
    name = name, outcome = outcome, model_type = "plain",
    train_raw = train_raw, prep_rec = prep_rec,
    x_cols = x_cols, final_fit = final_fit
  )
}

build_umap_hdb_bundle <- function(name, outcome, train_raw, prep_rec, final_fit,
                                  umap_model_file, hdb_fit, z_train,
                                  umap_settings) {
  x_cols <- bake(prep_rec, new_data = train_raw) %>%
    select(-any_of(c(outcome, "peiid"))) %>%
    mutate(across(where(is.logical), as.double)) %>%
    select(where(is.numeric)) %>%
    names()
  
  list(
    name = name, outcome = outcome, model_type = "umap_hdb",
    train_raw = train_raw, prep_rec = prep_rec, x_cols = x_cols,
    final_fit = final_fit, umap_model_file = umap_model_file,
    hdb_fit = hdb_fit, train_umap_features = as.matrix(z_train),
    umap_colnames = colnames(z_train), umap_settings = umap_settings
  )
}

# -----------------------------------------------------------------------------
# Evaluation: standard metrics + top-3 + hierarchical-parent accuracy
# -----------------------------------------------------------------------------
top_k_accuracy <- function(pred_tbl, truth_col, k = 3) {
  prob_cols <- names(pred_tbl)[
    str_detect(names(pred_tbl), "^\\.pred_") &
      names(pred_tbl) != ".pred_class"
  ]
  if (length(prob_cols) < k) return(NA_real_)
  
  prob_mat <- as.matrix(pred_tbl[prob_cols])
  class_names <- str_remove(prob_cols, "^\\.pred_")
  truth <- as.character(pred_tbl[[truth_col]])
  
  in_top_k <- vapply(seq_len(nrow(prob_mat)), function(i) {
    top_classes <- class_names[order(prob_mat[i, ], decreasing = TRUE)[1:k]]
    truth[i] %in% top_classes
  }, logical(1))
  
  mean(in_top_k, na.rm = TRUE)
}

evaluate_predictions <- function(pred_tbl, truth_col, k_for_topk = 3) {
  pr <- pred_tbl %>% rename(.truth = all_of(truth_col))
  out <- tibble(
    accuracy     = accuracy(pr,     truth = .truth, estimate = .pred_class)$.estimate,
    kap          = kap(pr,          truth = .truth, estimate = .pred_class)$.estimate,
    bal_accuracy = bal_accuracy(pr, truth = .truth, estimate = .pred_class)$.estimate,
    f_meas_macro = f_meas(pr,       truth = .truth, estimate = .pred_class,
                          estimator = "macro")$.estimate,
    top_k        = top_k_accuracy(pred_tbl, truth_col, k = k_for_topk),
    k            = k_for_topk
  )
  out
}

# -----------------------------------------------------------------------------
# Domain-shift diagnostic: train binary classifier on baked features
# -----------------------------------------------------------------------------
domain_shift_check <- function(ncss_baked, perry_baked, seed = 42) {
  common <- intersect(names(ncss_baked), names(perry_baked))
  ncss_b  <- ncss_baked  %>% select(all_of(common)) %>% mutate(.src = "ncss")
  perry_b <- perry_baked %>% select(all_of(common)) %>% mutate(.src = "perry")
  
  combined <- bind_rows(ncss_b, perry_b) %>%
    mutate(.src = factor(.src),
           across(where(is.numeric), ~ ifelse(is.finite(.x), .x, 0)))
  
  set.seed(seed)
  split <- initial_split(combined, prop = 0.8, strata = .src)
  tr <- training(split); te <- testing(split)
  
  rf_spec <- rand_forest(trees = 500, mtry = floor(sqrt(length(common))),
                         min_n = 5) %>%
    set_engine("ranger", probability = TRUE) %>%
    set_mode("classification")
  
  fit <- workflow() %>%
    add_model(rf_spec) %>%
    add_formula(.src ~ .) %>%
    fit(tr)
  
  pr <- predict(fit, te, type = "prob") %>%
    bind_cols(predict(fit, te, type = "class")) %>%
    bind_cols(te %>% select(.src))
  
  acc <- accuracy(pr, truth = .src, estimate = .pred_class)$.estimate
  
  imp <- workflows::extract_fit_engine(fit)$variable.importance
  if (is.null(imp)) {
    imp_tbl <- tibble(feature = character(), importance = numeric())
  } else {
    # For mtry classifier we set importance = "impurity" via default
    imp_tbl <- tibble(feature = names(imp), importance = as.numeric(imp)) %>%
      arrange(desc(importance))
  }
  
  list(
    accuracy = acc,
    n_features = length(common),
    top_features = imp_tbl %>% slice_head(n = 15),
    interpretation = case_when(
      acc < 0.6 ~ "Low domain shift: NCSS and Perry features are largely indistinguishable.",
      acc < 0.75 ~ "Moderate domain shift: features differ but not dramatically.",
      TRUE ~ "Strong domain shift: NCSS and Perry features are clearly different. Predictions should be interpreted with caution."
    )
  )
}

# -----------------------------------------------------------------------------
# Hierarchical-parent accuracy mapping (NCSS-style hierarchy)
# -----------------------------------------------------------------------------
HIERARCHY_PARENT <- list(
  taxonname_clean = "taxsubgrp_mod1",
  taxsubgrp_mod1  = "taxgrtgroup",
  taxgrtgroup     = "taxsuborder",
  taxsuborder     = "taxorder",
  taxorder        = NA_character_,
  taxpartsize     = NA_character_
)

# -----------------------------------------------------------------------------
# Perry deployment (carried over from your old script, lightly cleaned)
# -----------------------------------------------------------------------------
add_top1_margin <- function(pred_df) {
  pred_cols <- names(pred_df)[
    str_detect(names(pred_df), "^\\.pred_") & names(pred_df) != ".pred_class"
  ]
  pred_df %>%
    rowwise() %>%
    mutate(
      .probs  = list(c_across(all_of(pred_cols))),
      top1    = max(.probs, na.rm = TRUE),
      top2    = ifelse(length(.probs) > 1, sort(.probs, decreasing = TRUE)[2], NA_real_),
      margin  = ifelse(length(.probs) > 1, top1 - top2, NA_real_)
    ) %>%
    ungroup() %>%
    select(-.probs)
}

pad_like_train_raw <- function(new_df, train_df, outcome) {
  req <- setdiff(names(train_df), outcome)
  out <- new_df
  missing_cols <- setdiff(req, names(out))
  pad_log <- tibble(column = character(), pad_value = character())
  
  for (nm in missing_cols) {
    tmpl <- train_df[[nm]]
    if (str_detect(nm, "^prop_hzbase_")) {
      out[[nm]] <- 0; pad_log <- add_row(pad_log, column = nm, pad_value = "0")
    } else if (str_detect(nm, "^(any_|has_)") || is.logical(tmpl)) {
      out[[nm]] <- FALSE; pad_log <- add_row(pad_log, column = nm, pad_value = "FALSE")
    } else if (is.integer(tmpl)) {
      out[[nm]] <- NA_integer_; pad_log <- add_row(pad_log, column = nm, pad_value = "NA_int")
    } else if (is.numeric(tmpl)) {
      out[[nm]] <- NA_real_; pad_log <- add_row(pad_log, column = nm, pad_value = "NA_num")
    } else {
      out[[nm]] <- NA_character_; pad_log <- add_row(pad_log, column = nm, pad_value = "NA_chr")
    }
  }
  
  # Force every padded column to match the train template's type.
  # This catches the case where a column ALREADY exists in new_df but with
  # the wrong type (e.g., all-NA -> logical) -- the bug that broke POINT_X.
  for (nm in intersect(req, names(out))) {
    tmpl <- train_df[[nm]]
    if      (is.integer(tmpl))   out[[nm]] <- as.integer(out[[nm]])
    else if (is.logical(tmpl))   out[[nm]] <- as.logical(out[[nm]])
    else if (is.numeric(tmpl))   out[[nm]] <- as.numeric(out[[nm]])
    else if (is.character(tmpl)) out[[nm]] <- as.character(out[[nm]])
  }
  
  out <- out %>%
    mutate(across(where(is.logical), ~ tidyr::replace_na(.x, FALSE))) %>%
    select(any_of(c("peiid", req)))
  
  attr(out, "pad_log") <- pad_log
  out
}

deploy_plain_to_new <- function(bundle, new_pack) {
  outcome <- bundle$outcome
  new_raw <- pad_like_train_raw(new_pack, bundle$train_raw, outcome)
  
  common_cols <- intersect(names(new_raw), names(bundle$train_raw))
  for (nm in common_cols) {
    tmpl <- bundle$train_raw[[nm]]
    if      (is.integer(tmpl)) new_raw[[nm]] <- as.integer(new_raw[[nm]])
    else if (is.logical(tmpl)) new_raw[[nm]] <- as.logical(new_raw[[nm]])
    else if (is.numeric(tmpl)) new_raw[[nm]] <- as.numeric(new_raw[[nm]])
    else if (is.character(tmpl)) new_raw[[nm]] <- as.character(new_raw[[nm]])
    else if (is.factor(tmpl))
      new_raw[[nm]] <- factor(as.character(new_raw[[nm]]), levels = levels(tmpl))
  }
  
  bind_cols(
    new_raw %>% select(peiid),
    predict(bundle$final_fit, new_data = new_raw, type = "prob"),
    predict(bundle$final_fit, new_data = new_raw, type = "class")
  ) %>%
    add_top1_margin() %>%
    mutate(model_name = bundle$name, outcome = bundle$outcome,
           model_type = bundle$model_type)
}

deploy_umap_hdb_to_new <- function(bundle, new_pack) {
  outcome <- bundle$outcome
  new_raw <- pad_like_train_raw(new_pack, bundle$train_raw, outcome)
  
  x_new <- bake(bundle$prep_rec, new_data = new_raw) %>%
    select(-any_of(c(outcome, "peiid"))) %>%
    mutate(across(where(is.logical), as.double)) %>%
    select(where(is.numeric)) %>%
    select(all_of(bundle$x_cols)) %>%
    mutate(across(everything(), ~ ifelse(is.finite(.x), .x, 0)))
  
  z_new <- as.data.frame(uwot::umap_transform(as.matrix(x_new), bundle$umap_fit))
  names(z_new) <- bundle$umap_colnames
  
  train_z <- as.matrix(bundle$train_umap_features)
  train_cluster <- bundle$hdb_fit$cluster
  
  nearest_idx <- apply(as.matrix(z_new), 1, function(z) {
    d <- rowSums((train_z - matrix(z, nrow(train_z), ncol(train_z), byrow = TRUE))^2)
    which.min(d)
  })
  new_cluster <- train_cluster[nearest_idx]
  hdb_levels <- sort(unique(train_cluster))
  
  new_dat <- bind_cols(
    z_new,
    hdb_cluster = factor(new_cluster, levels = hdb_levels),
    hdb_noise   = factor(if_else(new_cluster == 0, "noise", "cluster"),
                         levels = c("cluster", "noise"))
  )
  
  bind_cols(
    new_raw %>% select(peiid),
    predict(bundle$final_fit, new_data = new_dat, type = "prob"),
    predict(bundle$final_fit, new_data = new_dat, type = "class")
  ) %>%
    add_top1_margin() %>%
    mutate(model_name = bundle$name, outcome = bundle$outcome,
           model_type = bundle$model_type)
}

deploy_bundle_to_new <- function(bundle, new_pack) {
  if (bundle$model_type == "plain") deploy_plain_to_new(bundle, new_pack)
  else if (bundle$model_type == "umap_hdb") deploy_umap_hdb_to_new(bundle, new_pack)
  else stop("Unknown model_type: ", bundle$model_type)
}