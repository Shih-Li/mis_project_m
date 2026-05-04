# ==============================================================================
# File: /script/04_compare_robust.R
# Purpose: Orchestrates the cross-comparison of parameter recovery between
#          OLS, Classical Diagnostics, MIS, LTS, and MM-estimators across the
#          full DGP grid (3 X distributions x 8 error distributions x 4
#          contamination topologies = 96 scenarios).
# Outputs: ../output/04_robust_comparison_results.rds
# ==============================================================================

# 1. Load Dependencies
library(dplyr)
library(purrr)
library(future)
library(furrr)
library(robustbase)

source("../R/dgp_factory.R")
source("../R/influence_injector.R")
source("../R/diagnostics_classical.R")
source("../R/fast_sens_topk.R")
source("../R/estimators_robust.R")
source("../R/sim_robust_engine.R")
source("../R/utils_checkpoint.R")

# 2. Global Configuration
sim_params <- list(
  n_iters   = 100,
  n_obs     = 1000,
  set_size  = 20,       # Contamination size (2% of data)
  magnitude = 10,       # Severe shift
  seed      = 20260503
)

set.seed(sim_params$seed)

# Setup Parallel Workers
num_workers <- max(1, parallel::detectCores() - 2)
cat(sprintf("Local environment, using %d workers.\n", num_workers))
plan(multisession, workers = num_workers)

# 3. Define the Full Reality Grid
#    3 X distributions x 8 error distributions x 4 outlier topologies = 96 cells
param_grid <- expand.grid(
  x_type         = c("normal", "mixed_normal", "contaminated"),
  error_type     = c("normal", "mixed_normal", "skewed_t", "golm",
                     "beta_logistic", "gpd", "contaminated", "pareto"),
  outlier_method = c("none", "vertical_outlier", "good_leverage", "bad_leverage"),
  stringsAsFactors = FALSE
)

cat(sprintf(
  paste0("Starting 04 Robust Comparison Suite.\n",
         "  X distributions:     %d\n",
         "  Error distributions: %d\n",
         "  Outlier topologies:  %d\n",
         "  Total scenarios:     %d\n",
         "  Iterations each:     %d\n",
         "  Total MC draws:      %d\n\n"),
  length(unique(param_grid$x_type)),
  length(unique(param_grid$error_type)),
  length(unique(param_grid$outlier_method)),
  nrow(param_grid),
  sim_params$n_iters,
  nrow(param_grid) * sim_params$n_iters
))

# 4. The Orchestrator Loop
for (i in seq_len(nrow(param_grid))) {
  
  p_current <- param_grid[i, , drop = FALSE]
  chunk_file <- sprintf("../output/temp_04/04_chunk_%03d.rds", i)
  
  if (is_computed(chunk_file)) {
    cat(sprintf("[%03d/%03d] Skipping (cached): x=%s | err=%s | out=%s\n",
                i, nrow(param_grid),
                p_current$x_type, p_current$error_type, p_current$outlier_method))
    next
  }
  
  cat(sprintf("[%03d/%03d] Computing: x=%s | err=%s | out=%s ... ",
              i, nrow(param_grid),
              p_current$x_type, p_current$error_type, p_current$outlier_method))
  
  # Run iterations in parallel
  scenario_results <- furrr::future_map_dfr(
    seq_len(sim_params$n_iters), function(iter_id) {
      
      tryCatch({
        run_robust_comparison_iter(
          iter           = iter_id,
          n              = sim_params$n_obs,
          p              = 1,
          x_type         = p_current$x_type,
          error_type     = p_current$error_type,
          outlier_method = p_current$outlier_method,
          k              = sim_params$set_size,
          magnitude      = sim_params$magnitude
        )
      }, error = function(e) {
        warning(sprintf("Iter %d failed for x=%s|err=%s|out=%s: %s",
                        iter_id, p_current$x_type, p_current$error_type,
                        p_current$outlier_method, e$message))
        return(NULL)
      })
      
    }, .options = furrr_options(seed = TRUE))
  
  # Safely save the chunk
  safe_save_rds(scenario_results, chunk_file)
  cat("Done.\n")
}

# 5. Final Assembly
cat("\nAll scenarios completed. Assembling final dataset...\n")
final_dataset <- compile_checkpoints(
  temp_dir          = "../output/temp_04",
  pattern           = "^04_chunk_.*\\.rds$",
  final_output_path = "../output/04_robust_comparison_results.rds",
  clear_temp        = FALSE
)

cat("\nScript 04 execution finished successfully.\n")

# ==============================================================================
# 6. Quick Sanity Check
# ==============================================================================
library(tidyr)

results <- readRDS("../output/04_robust_comparison_results.rds")

# 6a. Mean Absolute Bias
bias_table <- results %>%
  group_by(x_type, error_type, outlier_method) %>%
  summarise(across(starts_with("bias_"),
                   ~ mean(., na.rm = TRUE),
                   .names = "{.col}"),
            .groups = "drop") %>%
  rename_with(~ gsub("bias_", "", .x), starts_with("bias_"))

cat("\n--- Mean Absolute Bias ---\n")
print(bias_table, n = Inf)

# 6b. 95% CI Coverage
coverage_table <- results %>%
  group_by(x_type, error_type, outlier_method) %>%
  summarise(across(starts_with("cov_"),
                   ~ mean(., na.rm = TRUE) * 100,
                   .names = "{.col}"),
            .groups = "drop") %>%
  rename_with(~ gsub("cov_", "", .x), starts_with("cov_"))

cat("\n--- 95% CI Coverage (%) ---\n")
print(coverage_table, n = Inf)

# 6c. RMSE
rmse_table <- results %>%
  group_by(x_type, error_type, outlier_method) %>%
  summarise(across(starts_with("coef_"),
                   ~ sqrt(mean((. - 1)^2, na.rm = TRUE)),
                   .names = "{.col}"),
            .groups = "drop") %>%
  rename_with(~ gsub("coef_", "", .x), starts_with("coef_"))

cat("\n--- Root Mean Square Error (RMSE) ---\n")
print(rmse_table, n = Inf)