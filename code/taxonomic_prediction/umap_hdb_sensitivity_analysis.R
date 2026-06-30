# ---
#   title: "Test_UMAP_tune_05_14_2026"
# format: html
# ---
#   

# =============================================================================
# umap_hdb_sensitivity_analysis.R
#
# Tests whether UMAP+HDBSCAN underperformance is due to fixed default hyperparams
# or whether it's a structural property of the representation.
#
# For ONE outcome:
#   1. Build a small UMAP × HDBSCAN parameter grid
#   2. For each (n_neighbors, min_dist, minPts) combination:
#       - Fit UMAP+HDB on canonical 80/20 NCSS train split
#       - Tune RF and XGB on the UMAP representation (inner CV)
#       - Score on held-out NCSS test split
#       - Deploy to Perry, score there
#   3. Compare to canonical default-UMAP baseline
#
# Output: CSV with kappa per (outcome × engine × UMAP params) for both
# NCSS test and Perry external, plus tuning grid response surface plots.
# =============================================================================

library(tidyverse)
library(here)
library(yardstick)

source(here("code", "modeling_functions_05_02_2026.R"))
source(here("code", "pedon_pack_functions_05_02_2026.R"))

# -----------------------------------------------------------------------------
# CONFIG — adjust these as needed
# -----------------------------------------------------------------------------
OUTCOMES_TO_TEST <- c("taxorder", "taxsuborder", "taxpartsize", "taxonname_clean")  
GRID_SIZE_INNER <- 6        # smaller inner grid to keep time tractable
INNER_CV_V       <- 3        # 3 inner folds instead of 5
SEED             <- 2026

# UMAP × HDBSCAN parameter grid
umap_grid <- expand_grid(
  n_neighbors  = c(5, 15, 30, 50),
  min_dist     = c(0.01, 0.1, 0.3),
  minPts       = c(5, 10, 20)
  # n_components fixed at 10 (default)
)
cat("UMAP grid size:", nrow(umap_grid), "combinations per outcome × engine\n")
cat("Total combinations:", nrow(umap_grid) * length(OUTCOMES_TO_TEST) * 2,
    "(2 engines)\n\n")

MIN_N_MAP <- c(taxorder = 15, taxsuborder = 15, taxgrtgroup = 15,
               taxsubgrp_mod1 = 12, taxpartsize = 20,
               taxonname_clean = 5)

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------
ncss_pack_files <- list.files(here("outputs"),
                              pattern = "^pedon_pack_150_.*\\.csv$",
                              full.names = TRUE)
ncss_pack_path <- ncss_pack_files[which.max(file.mtime(ncss_pack_files))]
ncss_full <- read.csv(ncss_pack_path) %>%
  mutate(peiid = as.character(peiid))
cat("Loaded NCSS pack:", basename(ncss_pack_path), "—", nrow(ncss_full), "pedons\n")

perry_pack_files <- list.files(here("outputs"),
                               pattern = "^perry_pack_150_.*\\.csv$",
                               full.names = TRUE)
perry_pack_path <- perry_pack_files[which.max(file.mtime(perry_pack_files))]
perry_pack <- read.csv(perry_pack_path) %>%
  mutate(peiid = as.character(peiid),
         POINT_X = as.numeric(POINT_X),
         POINT_Y = as.numeric(POINT_Y))
cat("Loaded Perry pack:", basename(perry_pack_path), "—", nrow(perry_pack), "pedons\n")

perry_fp_files <- list.files(here("outputs"),
                             pattern = "^perry_fp_clean_.*\\.csv$",
                             full.names = TRUE)
perry_fp_path <- perry_fp_files[which.max(file.mtime(perry_fp_files))]
perry_fp <- read.csv(perry_fp_path) %>%
  mutate(peiid = as.character(Point))

# Perry truth-builder (same as in your main pipeline)
build_perry_truth <- function(perry_fp_clean, truth_col) {
  if (!truth_col %in% names(perry_fp_clean)) return(tibble())
  dat <- perry_fp_clean %>%
    mutate(truth_val = na_if(str_squish(as.character(.data[[truth_col]])), ""))
  if (truth_col == "taxpartsize") {
    dat <- dat %>% mutate(truth_val = str_to_lower(truth_val))
  }
  dat %>%
    group_by(peiid) %>%
    summarise(
      n_truth_unique = n_distinct(truth_val[!is.na(truth_val)]),
      truth = {
        nv <- truth_val[!is.na(truth_val)]
        if (length(nv) == 0) NA_character_ else
          names(sort(table(nv), decreasing = TRUE))[1]
      },
      .groups = "drop"
    ) %>%
    filter(!is.na(truth), n_truth_unique == 1)
}

# -----------------------------------------------------------------------------
# CORE FUNCTION: fit one UMAP+HDB configuration and score on NCSS test + Perry
# -----------------------------------------------------------------------------
fit_and_score_one_umap_config <- function(outcome, engine,
                                          n_neighbors, min_dist, minPts,
                                          ncss_train, ncss_test, perry_pack,
                                          perry_truth) {
  
  cat(sprintf("  outcome=%s engine=%s nn=%d md=%.2f mpts=%d ... ",
              outcome, engine, n_neighbors, min_dist, minPts))
  start_time <- Sys.time()
  
  umap_settings <- list(
    n_components = 10,
    n_neighbors  = n_neighbors,
    min_dist     = min_dist,
    metric       = "euclidean",
    minPts       = minPts
  )
  
  # Try to fit; if it crashes, return NA row
  result <- tryCatch({
    res <- fit_one_model(
      train_df = ncss_train,
      test_df  = ncss_test,
      outcome  = outcome,
      engine   = engine,
      space    = "umap_hdb",
      inner_v  = INNER_CV_V,
      grid_size = GRID_SIZE_INNER,
      seed     = SEED,
      umap_settings = umap_settings
    )
    
    # Score on NCSS test
    pr_ncss <- res$pred %>% rename(.truth = all_of(outcome))
    ncss_kap <- kap(pr_ncss, truth = .truth, estimate = .pred_class)$.estimate
    ncss_acc <- accuracy(pr_ncss, truth = .truth, estimate = .pred_class)$.estimate
    
    # Build a bundle for Perry deployment
    # Note: we need umap_model_file. Save the umap fit to a temp file
    tmp_umap <- tempfile(pattern = paste0("umap_sens_", outcome, "_"),
                         fileext = ".uwot")
    uwot::save_uwot(res$umap_fit, tmp_umap)
    
    bundle <- build_umap_hdb_bundle(
      name = paste0(outcome, "_", engine, "_umap_nn", n_neighbors,
                    "_md", min_dist, "_mpts", minPts),
      outcome = outcome,
      train_raw = res$train_raw,
      prep_rec = res$prep_rec,
      final_fit = res$final_fit,
      umap_model_file = tmp_umap,
      hdb_fit = res$hdb_fit,
      z_train = res$z_train,
      umap_settings = umap_settings
    )
    bundle$umap_fit <- res$umap_fit  # attach for deployment
    
    # Deploy to Perry
    perry_pred <- tryCatch(
      deploy_bundle_to_new(bundle, perry_pack) %>%
        mutate(peiid = as.character(peiid)),
      error = function(e) {
        cat("PERRY DEPLOY ERROR: ", e$message, " ")
        NULL
      }
    )
    
    perry_kap <- NA_real_
    perry_acc <- NA_real_
    perry_n   <- 0
    
    if (!is.null(perry_pred)) {
      joined <- perry_pred %>%
        left_join(perry_truth %>% select(peiid, truth),
                  by = "peiid") %>%
        filter(!is.na(truth))
      
      if (nrow(joined) >= 2) {
        # Build harmonized factor levels
        lvl <- union(
          as.character(unique(joined$truth)),
          as.character(unique(joined$.pred_class))
        )
        joined <- joined %>%
          mutate(
            truth_fct = factor(truth, levels = lvl),
            pred_fct  = factor(as.character(.pred_class), levels = lvl)
          )
        perry_kap <- kap(joined, truth = truth_fct, estimate = pred_fct)$.estimate
        perry_acc <- accuracy(joined, truth = truth_fct, estimate = pred_fct)$.estimate
        perry_n <- nrow(joined)
      }
    }
    
    unlink(tmp_umap)  # clean up temp file
    
    tibble(
      outcome = outcome, engine = engine,
      n_neighbors = n_neighbors, min_dist = min_dist, minPts = minPts,
      ncss_kap = ncss_kap, ncss_acc = ncss_acc,
      perry_kap = perry_kap, perry_acc = perry_acc, perry_n = perry_n,
      status = "ok"
    )
  }, error = function(e) {
    cat("FIT ERROR: ", e$message, " ")
    tibble(
      outcome = outcome, engine = engine,
      n_neighbors = n_neighbors, min_dist = min_dist, minPts = minPts,
      ncss_kap = NA_real_, ncss_acc = NA_real_,
      perry_kap = NA_real_, perry_acc = NA_real_, perry_n = 0L,
      status = paste("error:", e$message)
    )
  })
  
  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
  cat(sprintf("done [%.0fs] ncss_kap=%.3f perry_kap=%.3f\n",
              elapsed, result$ncss_kap, result$perry_kap))
  
  result
}

# -----------------------------------------------------------------------------
# RUN THE GRID FOR EACH OUTCOME
# -----------------------------------------------------------------------------
all_sensitivity_results <- list()

for (outcome in OUTCOMES_TO_TEST) {
  cat("\n=================================================\n")
  cat("OUTCOME:", outcome, "\n")
  cat("=================================================\n")
  
  # Prepare NCSS data: drop other taxa, filter by min_n, create canonical 80/20 split
  min_n <- MIN_N_MAP[[outcome]]
  df_prepped <- ncss_full %>%
    select(-any_of(setdiff(TAXONOMY_DROP_COLS, outcome))) %>%
    prep_outcome(outcome, min_n = min_n)
  
  set.seed(2026)  # canonical seed (no rep multiplier — this is single-split)
  s <- initial_split(df_prepped, prop = 0.8, strata = !!sym(outcome))
  ncss_train <- training(s)
  ncss_test  <- testing(s)
  
  cat("NCSS train:", nrow(ncss_train), " test:", nrow(ncss_test),
      " classes:", nlevels(ncss_train[[outcome]]), "\n\n")
  
  # Perry truth for this outcome
  truth_col <- if (outcome == "taxonname_clean") "series_clean" else outcome
  perry_truth <- build_perry_truth(perry_fp, truth_col)
  cat("Perry truth pedons:", nrow(perry_truth), "\n\n")
  
  # Loop over grid × engines
  for (engine in c("rf", "xgb")) {
    cat(sprintf("\n--- ENGINE: %s ---\n", engine))
    for (i in seq_len(nrow(umap_grid))) {
      row <- umap_grid[i, ]
      result <- fit_and_score_one_umap_config(
        outcome = outcome,
        engine = engine,
        n_neighbors = row$n_neighbors,
        min_dist = row$min_dist,
        minPts = row$minPts,
        ncss_train = ncss_train,
        ncss_test = ncss_test,
        perry_pack = perry_pack,
        perry_truth = perry_truth
      )
      all_sensitivity_results <- append(all_sensitivity_results, list(result))
    }
  }
}

# -----------------------------------------------------------------------------
# COMBINE AND SAVE
# -----------------------------------------------------------------------------
sensitivity_df <- bind_rows(all_sensitivity_results)

cat("\n\n=================================================\n")
cat("SENSITIVITY ANALYSIS COMPLETE\n")
cat("=================================================\n")
print(sensitivity_df, n = Inf)

write.csv(sensitivity_df,
          here("outputs", "reports_05_02", "umap_sensitivity_results.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# SUMMARY: BEST PER OUTCOME × ENGINE × DATASET
# -----------------------------------------------------------------------------
cat("\n\n=================================================\n")
cat("BEST CONFIG PER OUTCOME × ENGINE (by NCSS kap)\n")
cat("=================================================\n")
best_by_ncss <- sensitivity_df %>%
  filter(status == "ok") %>%
  group_by(outcome, engine) %>%
  slice_max(ncss_kap, n = 1, with_ties = FALSE) %>%
  ungroup()
print(best_by_ncss)

cat("\n\n=================================================\n")
cat("BEST CONFIG PER OUTCOME × ENGINE (by Perry kap)\n")
cat("=================================================\n")
best_by_perry <- sensitivity_df %>%
  filter(status == "ok") %>%
  group_by(outcome, engine) %>%
  slice_max(perry_kap, n = 1, with_ties = FALSE) %>%
  ungroup()
print(best_by_perry)

# -----------------------------------------------------------------------------
# COMPARE TO YOUR EXISTING DEFAULT-UMAP BASELINE
# -----------------------------------------------------------------------------
cat("\n\n=================================================\n")
cat("COMPARISON TO DEFAULTS (n_neighbors=15, min_dist=0.1, minPts=10)\n")
cat("=================================================\n")
defaults <- sensitivity_df %>%
  filter(status == "ok",
         n_neighbors == 15, min_dist == 0.1, minPts == 10)

cat("\nDefault-config results:\n")
print(defaults)

cat("\nBest-found - Default kappa improvements:\n")
improvement <- best_by_ncss %>%
  select(outcome, engine, best_ncss_kap = ncss_kap,
         best_perry_kap = perry_kap) %>%
  left_join(
    defaults %>% select(outcome, engine,
                        default_ncss_kap = ncss_kap,
                        default_perry_kap = perry_kap),
    by = c("outcome", "engine")
  ) %>%
  mutate(
    ncss_improvement = best_ncss_kap - default_ncss_kap,
    perry_improvement = best_perry_kap - default_perry_kap
  )
print(improvement)

# Save the improvement table — this is what goes in the methods/results
write.csv(improvement,
          here("outputs", "reports_05_02", "umap_sensitivity_improvement_05_14_26.csv"),
          row.names = FALSE)

cat("\nSensitivity analysis files saved to outputs/reports_05_02/\n")
# Rename for clarity
file.rename(
  here("outputs", "reports_05_02", "umap_sensitivity_results.csv"),
  here("outputs", "reports_05_02", "umap_sensitivity_results_part2_others.csv")
)
library(tidyverse)
library(here)

part1 <- read.csv(here("outputs", "reports_05_02",
                       "umap_sensitivity_results_part1_gg_subgrp.csv"))
part2 <- read.csv(here("outputs", "reports_05_02",
                       "umap_sensitivity_results_part2_others.csv"))

sensitivity_full <- bind_rows(part1, part2)

write.csv(sensitivity_full,
          here("outputs", "reports_05_02", "umap_sensitivity_results_FULL.csv"),
          row.names = FALSE)

cat("Combined results: ", nrow(sensitivity_full), "rows across",
    n_distinct(sensitivity_full$outcome), "outcomes\n")
best_umap_settings <- sensitivity_full %>%
  filter(status == "ok") %>%
  group_by(outcome, engine) %>%
  slice_max(ncss_kap, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(outcome, engine, n_neighbors, min_dist, minPts,
         ncss_kap_best = ncss_kap, perry_kap_best = perry_kap)

write.csv(best_umap_settings,
          here("outputs", "reports_05_02", "best_umap_settings.csv"),
          row.names = FALSE)

print(best_umap_settings, n = Inf)

# Load existing NCSS + Perry results
existing_results <- read.csv(here("outputs", "reports_05_02",
                                  "combined_ncss_perry_results.csv"))

# Build "tuned UMAP" rows from sensitivity analysis
tuned_umap_rows <- sensitivity_full %>%
  filter(status == "ok") %>%
  group_by(outcome, engine) %>%
  slice_max(ncss_kap, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # Two rows per outcome × engine — one for NCSS, one for Perry
  pivot_longer(cols = c(ncss_kap, perry_kap, ncss_acc, perry_acc),
               names_to = c("dataset_short", "metric"),
               names_sep = "_") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate(
    dataset = if_else(dataset_short == "ncss",
                      "NCSS (internal test)", "Perry (external)"),
    space = "umap_hdb_tuned",   # distinct from plain umap_hdb
    is_min5 = FALSE,
    bundle = paste0(outcome, "_", engine, "_umap_tuned_nn",
                    n_neighbors, "_md", min_dist, "_mpts", minPts),
    bal_accuracy = NA_real_,
    f_meas_macro = NA_real_,
    top_k = NA_real_,
    k = NA_real_,
    n_pedons = NA_integer_
  ) %>%
  select(dataset, outcome, engine, space, is_min5, bundle,
         accuracy = acc, kap, bal_accuracy, f_meas_macro,
         top_k, k, n_pedons)

# Combine
combined_with_tuned <- bind_rows(existing_results, tuned_umap_rows)

write.csv(combined_with_tuned,
          here("outputs", "reports_05_02",
               "combined_ncss_perry_with_tuned_umap.csv"),
          row.names = FALSE)
