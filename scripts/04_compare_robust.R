# ==============================================================================
# File: /scripts/04_compare_robust.R
# Purpose: Orchestrate the robust-estimator comparison using exact Dinkelbach
#          MIS detection and Selection-Adjusted Permutation (SAP) calibration.
#
# Estimators (10):
#   full, cd, lev, dfb, mis_alpha, mis_oracle, mis_peel, mis_sap, mm, lts
#
# Outputs:
#   ../output/04sap_robust_comparison_results.rds
#   ../output/04sap_summary_tables.rds
#   ../output/04sap_bias_distributional.rds
#
# Run with the working directory set to /scripts.
# ==============================================================================


# ==============================================================================
# 1. Packages and project functions
# ==============================================================================

library(dplyr)
library(future)
library(furrr)
library(robustbase)

# helpers_local.R must be sourced first.
source("../R/helpers_local.R")
source("../R/dgp_factory.R")
source("../R/influence_injector.R")
source("../R/diagnostics_classical.R")
source("../R/estimators_robust.R")
source("../R/dynamic_k_adaptive.R")
source("../R/dinkelbach_topk.R")
source("../R/leverage_k.R")
source("../R/iterative_peel_v2.R")
source("../R/iterative_peel_sap.R")
source("../R/sim_robust_engine_v2.R")
source("../R/utils_checkpoint.R")

# Guard against accidentally using an old engine left in memory.
required_engine_args <- c(
  "sap_alpha", "sap_B_perm", "sap_k_grid", "sap_max_iter"
)
missing_engine_args <- setdiff(
  required_engine_args,
  names(formals(run_robust_comparison_iter_v2))
)
if (length(missing_engine_args) > 0L) {
  stop(
    "Old robust engine loaded. Missing SAP arguments: ",
    paste(missing_engine_args, collapse = ", "),
    ". Re-source ../R/sim_robust_engine_v2.R."
  )
}

required_functions <- c(
  "generate_complex_data", "apply_influence_shift", "get_classical_set",
  "fit_clean_ols", "fit_mm_estimator", "fit_lts_estimator",
  "alpha_k", "oracle_k", "dinkelbach_topk", "dinkelbach_topk_lm",
  "iterative_peel_v2", "iterative_peel_sap",
  "run_robust_comparison_iter_v2",
  "safe_save_rds", "is_computed", "compile_checkpoints"
)
missing_functions <- required_functions[
  !vapply(required_functions, exists, logical(1), mode = "function")
]
if (length(missing_functions) > 0L) {
  stop(
    "Missing project functions: ",
    paste(missing_functions, collapse = ", ")
  )
}


# ==============================================================================
# 2. Configuration
# ==============================================================================

sim_params <- list(
  # 1L = smoke test, 10L = pilot, 100L = final simulation.
  n_iters = 100L,
  n_obs = 5000L,
  set_size = 50L,
  magnitude = 10,
  seed = 20260503L,
  
  sap_alpha = 0.05,
  sap_B_perm = 199L,
  sap_k_grid = c(1L, 2L, 5L, 10L, 20L, 50L, 100L),
  sap_max_iter = 1L
)

param_grid <- expand.grid(
  x_type = c("normal", "mixed_normal", "contaminated"),
  error_type = c(
    "normal", "mixed_normal", "skewed_t", "golm",
    "beta_logistic", "gpd", "contaminated", "pareto"
  ),
  outlier_method = c(
    "none", "vertical_outlier", "good_leverage", "bad_leverage"
  ),
  stringsAsFactors = FALSE
)

set.seed(sim_params$seed)

n_scenarios <- nrow(param_grid)
n_expected <- n_scenarios * sim_params$n_iters
num_workers <- min(
  max(1L, future::availableCores() - 2L),
  sim_params$n_iters
)

cat(sprintf("Local environment: using %d worker(s).\n", num_workers))
future::plan(future::multisession, workers = num_workers)

cat(sprintf(
  paste0(
    "\nStarting Script 04: Robust MIS-SAP Comparison\n",
    "  Scenarios:           %d\n",
    "  Iterations/scenario: %d\n",
    "  Expected draws:      %d\n",
    "  Sample size:         %d\n",
    "  Injected k:          %d\n",
    "  SAP permutations:    %d\n",
    "  SAP k-grid:          %s\n\n"
  ),
  n_scenarios,
  sim_params$n_iters,
  n_expected,
  sim_params$n_obs,
  sim_params$set_size,
  sim_params$sap_B_perm,
  paste(sim_params$sap_k_grid, collapse = ", ")
))


# ==============================================================================
# 3. Checkpointed parallel simulation
# ==============================================================================

for (i in seq_len(n_scenarios)) {
  p_current <- param_grid[i, , drop = FALSE]
  chunk_file <- sprintf(
    "../output/temp_04sap/04sap_chunk_%03d.rds", i
  )
  
  if (is_computed(chunk_file)) {
    cat(sprintf(
      "[%03d/%03d] Cached: x=%s | error=%s | outlier=%s\n",
      i, n_scenarios,
      p_current$x_type,
      p_current$error_type,
      p_current$outlier_method
    ))
    next
  }
  
  cat(sprintf(
    "[%03d/%03d] Running: x=%s | error=%s | outlier=%s ... ",
    i, n_scenarios,
    p_current$x_type,
    p_current$error_type,
    p_current$outlier_method
  ))
  
  scenario_results <- furrr::future_map_dfr(
    seq_len(sim_params$n_iters),
    function(iter_id) {
      tryCatch(
        run_robust_comparison_iter_v2(
          iter = iter_id,
          n = sim_params$n_obs,
          p = 1L,
          x_type = p_current$x_type,
          error_type = p_current$error_type,
          outlier_method = p_current$outlier_method,
          k = sim_params$set_size,
          magnitude = sim_params$magnitude,
          sap_alpha = sim_params$sap_alpha,
          sap_B_perm = sim_params$sap_B_perm,
          sap_k_grid = sim_params$sap_k_grid,
          sap_max_iter = sim_params$sap_max_iter
        ),
        error = function(e) {
          warning(sprintf(
            "Iter %d failed for x=%s|error=%s|outlier=%s: %s",
            iter_id,
            p_current$x_type,
            p_current$error_type,
            p_current$outlier_method,
            conditionMessage(e)
          ))
          NULL
        }
      )
    },
    # Scenario-specific seeds remain stable when completed chunks are skipped.
    .options = furrr::furrr_options(seed = sim_params$seed + i)
  )
  
  safe_save_rds(scenario_results, chunk_file)
  
  cat(sprintf(
    "Done: %d successful, %d failed.\n",
    nrow(scenario_results),
    sim_params$n_iters - nrow(scenario_results)
  ))
}

future::plan(future::sequential)


# ==============================================================================
# 4. Final assembly and integrity checks
# ==============================================================================

cat("\nAssembling checkpoint files...\n")

compile_checkpoints(
  temp_dir = "../output/temp_04sap",
  pattern = "^04sap_chunk_.*\\.rds$",
  final_output_path = "../output/04sap_robust_comparison_results.rds",
  clear_temp = FALSE
)

results <- readRDS("../output/04sap_robust_comparison_results.rds")

cat(sprintf(
  "Observed rows: %d | Expected rows: %d | Missing rows: %d\n",
  nrow(results), n_expected, n_expected - nrow(results)
))

if (nrow(results) != n_expected) {
  warning(
    "The final dataset does not contain the expected number of rows. ",
    "Review warnings and checkpoint sizes before final analysis."
  )
}

duplicate_rows <- results %>%
  count(x_type, error_type, outlier_method, iter, name = "n") %>%
  filter(n > 1L)

if (nrow(duplicate_rows) > 0L) {
  warning("Duplicate scenario-iteration rows detected.")
  print(duplicate_rows, n = Inf)
}


# ==============================================================================
# 5. Health checks
# ==============================================================================

cat("\n--- Method Health Checks ---\n")

health_table <- results %>%
  summarise(
    n_rows = n(),
    sap_error_rate = mean(peel_sap_stop == "error", na.rm = TRUE),
    peel_v2_error_rate = mean(peel_v2_stop == "error", na.rm = TRUE),
    mm_convergence_rate = mean(mm_converged, na.rm = TRUE),
    sap_missing_p_rate = mean(!is.finite(peel_sap_final_p)),
    sap_missing_coef_rate = mean(!is.finite(coef_mis_sap))
  )

print(health_table)


# ==============================================================================
# 6. Selected-set sizes and stopping diagnostics
# ==============================================================================

cat("\n--- Selected or Flagged k ---\n")

k_table <- results %>%
  group_by(outlier_method) %>%
  summarise(
    true_k = mean(set_size, na.rm = TRUE),
    cd_mean = mean(k_cd, na.rm = TRUE),
    lev_mean = mean(k_lev, na.rm = TRUE),
    dfb_mean = mean(k_dfb, na.rm = TRUE),
    alpha_mean = mean(k_alpha, na.rm = TRUE),
    oracle_mean = mean(k_oracle, na.rm = TRUE),
    peel_v2_mean = mean(k_peel_v2, na.rm = TRUE),
    peel_v2_med = median(k_peel_v2, na.rm = TRUE),
    sap_mean = mean(k_peel_sap, na.rm = TRUE),
    sap_med = median(k_peel_sap, na.rm = TRUE),
    sap_q25 = unname(quantile(k_peel_sap, 0.25, na.rm = TRUE)),
    sap_q75 = unname(quantile(k_peel_sap, 0.75, na.rm = TRUE)),
    .groups = "drop"
  )

print(k_table)


cat("\n--- Peel-v2 Stop Reasons ---\n")

peel_v2_stops <- results %>%
  count(outlier_method, peel_v2_stop, name = "n") %>%
  group_by(outlier_method) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  arrange(outlier_method, desc(n))

print(peel_v2_stops, n = Inf)


cat("\n--- SAP Stop Reasons ---\n")

sap_stops <- results %>%
  count(outlier_method, peel_sap_stop, name = "n") %>%
  group_by(outlier_method) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  arrange(outlier_method, desc(n))

print(sap_stops, n = Inf)


cat("\n--- SAP Process Summary ---\n")

sap_diagnostics <- results %>%
  group_by(outlier_method) %>%
  summarise(
    detection_rate = mean(k_peel_sap > 0L, na.rm = TRUE),
    exact_k_rate = mean(k_peel_sap == set_size, na.rm = TRUE),
    k_mean = mean(k_peel_sap, na.rm = TRUE),
    k_med = median(k_peel_sap, na.rm = TRUE),
    iterations_mean = mean(peel_sap_iters, na.rm = TRUE),
    final_p_mean = mean(peel_sap_final_p, na.rm = TRUE),
    final_p_med = median(peel_sap_final_p, na.rm = TRUE),
    minimum_p_mean = mean(peel_sap_min_p, na.rm = TRUE),
    peak_excess_mean = mean(peel_sap_peak_excess, na.rm = TRUE),
    peak_excess_med = median(peel_sap_peak_excess, na.rm = TRUE),
    .groups = "drop"
  )

print(sap_diagnostics)


# ==============================================================================
# 7. Detection overlap
# ==============================================================================

cat("\n--- Detection Overlap: Mean (SD) ---\n")

overlap_table <- results %>%
  filter(outlier_method != "none") %>%
  group_by(outlier_method) %>%
  summarise(
    cd_mean = mean(overlap_cd, na.rm = TRUE),
    cd_sd = sd(overlap_cd, na.rm = TRUE),
    lev_mean = mean(overlap_lev, na.rm = TRUE),
    lev_sd = sd(overlap_lev, na.rm = TRUE),
    dfb_mean = mean(overlap_dfb, na.rm = TRUE),
    dfb_sd = sd(overlap_dfb, na.rm = TRUE),
    alpha_mean = mean(overlap_mis_alpha, na.rm = TRUE),
    alpha_sd = sd(overlap_mis_alpha, na.rm = TRUE),
    oracle_mean = mean(overlap_mis_oracle, na.rm = TRUE),
    oracle_sd = sd(overlap_mis_oracle, na.rm = TRUE),
    peel_v2_mean = mean(overlap_peel_v2, na.rm = TRUE),
    peel_v2_sd = sd(overlap_peel_v2, na.rm = TRUE),
    sap_mean = mean(overlap_peel_sap, na.rm = TRUE),
    sap_sd = sd(overlap_peel_sap, na.rm = TRUE),
    .groups = "drop"
  )

print(overlap_table)


# ==============================================================================
# 8. Coverage
# ==============================================================================

cat("\n--- 95% CI Coverage by Simulation Cell (%) ---\n")

coverage_table <- results %>%
  group_by(x_type, error_type, outlier_method) %>%
  summarise(
    n_iter = n(),
    across(
      starts_with("cov_"),
      ~ mean(.x, na.rm = TRUE) * 100
    ),
    .groups = "drop"
  )

print(coverage_table, n = Inf)


cat("\n--- Coverage by Contamination: Mean% (Monte Carlo SE%) ---\n")

coverage_summary <- results %>%
  group_by(outlier_method) %>%
  summarise(
    across(
      starts_with("cov_"),
      ~ {
        p_hat <- mean(.x, na.rm = TRUE)
        n_valid <- sum(!is.na(.x))
        if (n_valid == 0L) {
          NA_character_
        } else {
          mc_se <- sqrt(p_hat * (1 - p_hat) / n_valid)
          sprintf("%.1f (%.1f)", 100 * p_hat, 100 * mc_se)
        }
      }
    ),
    .groups = "drop"
  )

print(coverage_summary)


# ==============================================================================
# 9. Absolute bias and RMSE
# ==============================================================================

cat("\n--- Absolute Bias by Simulation Cell ---\n")

bias_distributional <- results %>%
  group_by(x_type, error_type, outlier_method) %>%
  summarise(
    across(
      starts_with("bias_"),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),
        med = ~ median(.x, na.rm = TRUE),
        q25 = ~ unname(quantile(.x, 0.25, na.rm = TRUE)),
        q75 = ~ unname(quantile(.x, 0.75, na.rm = TRUE))
      ),
      .names = "{.col}__{.fn}"
    ),
    .groups = "drop"
  )


cat("\n--- Absolute Bias by Contamination: Mean [Median] ---\n")

bias_summary <- results %>%
  group_by(outlier_method) %>%
  summarise(
    across(
      starts_with("bias_"),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        med = ~ median(.x, na.rm = TRUE)
      ),
      .names = "{.col}__{.fn}"
    ),
    .groups = "drop"
  )

print(bias_summary)


cat("\n--- RMSE by Simulation Cell ---\n")

rmse_table <- results %>%
  group_by(x_type, error_type, outlier_method) %>%
  summarise(
    across(
      starts_with("coef_"),
      ~ sqrt(mean((.x - 1)^2, na.rm = TRUE)),
      .names = "rmse_{.col}"
    ),
    .groups = "drop"
  )

print(rmse_table, n = Inf)


# ==============================================================================
# 10. SAP versus peel-v2
# ==============================================================================

cat("\n--- SAP versus Peel-v2 ---\n")

sap_vs_peel <- results %>%
  filter(outlier_method != "none") %>%
  group_by(outlier_method) %>%
  summarise(
    k_peel_v2_mean = mean(k_peel_v2, na.rm = TRUE),
    k_peel_v2_med = median(k_peel_v2, na.rm = TRUE),
    k_sap_mean = mean(k_peel_sap, na.rm = TRUE),
    k_sap_med = median(k_peel_sap, na.rm = TRUE),
    overlap_peel_v2 = mean(overlap_peel_v2, na.rm = TRUE),
    overlap_sap = mean(overlap_peel_sap, na.rm = TRUE),
    bias_peel_v2 = mean(bias_mis_peel, na.rm = TRUE),
    bias_sap = mean(bias_mis_sap, na.rm = TRUE),
    coverage_peel_v2 = mean(cov_mis_peel, na.rm = TRUE) * 100,
    coverage_sap = mean(cov_mis_sap, na.rm = TRUE) * 100,
    cpu_peel_v2_med = median(cpu_peel_v2, na.rm = TRUE),
    cpu_sap_med = median(cpu_peel_sap, na.rm = TRUE),
    .groups = "drop"
  )

print(sap_vs_peel)


# ==============================================================================
# 11. Complete method runtime
# ==============================================================================

cat("\n--- Complete Runtime: Median [Q25, Q75] Seconds ---\n")

runtime_table <- results %>%
  summarise(
    across(
      starts_with("cpu_"),
      list(
        med = ~ median(.x, na.rm = TRUE),
        q25 = ~ unname(quantile(.x, 0.25, na.rm = TRUE)),
        q75 = ~ unname(quantile(.x, 0.75, na.rm = TRUE))
      ),
      .names = "{.col}__{.fn}"
    )
  )

print(t(runtime_table))


# ==============================================================================
# 12. Key SAP result: bad leverage by error distribution
# ==============================================================================

cat("\n--- SAP under Bad Leverage by Error Distribution ---\n")

sap_bad_leverage <- results %>%
  filter(outlier_method == "bad_leverage") %>%
  group_by(error_type) %>%
  summarise(
    detection_rate = mean(k_peel_sap > 0L, na.rm = TRUE),
    exact_k_rate = mean(k_peel_sap == set_size, na.rm = TRUE),
    k_mean = mean(k_peel_sap, na.rm = TRUE),
    k_med = median(k_peel_sap, na.rm = TRUE),
    overlap_mean = mean(overlap_peel_sap, na.rm = TRUE),
    bias_mean = mean(bias_mis_sap, na.rm = TRUE),
    bias_med = median(bias_mis_sap, na.rm = TRUE),
    coverage_sap = mean(cov_mis_sap, na.rm = TRUE) * 100,
    coverage_peel_v2 = mean(cov_mis_peel, na.rm = TRUE) * 100,
    coverage_mm = mean(cov_mm, na.rm = TRUE) * 100,
    coverage_oracle = mean(cov_mis_oracle, na.rm = TRUE) * 100,
    final_p_med = median(peel_sap_final_p, na.rm = TRUE),
    peak_excess_med = median(peel_sap_peak_excess, na.rm = TRUE),
    .groups = "drop"
  )

print(sap_bad_leverage)


# ==============================================================================
# 13. Heavy-tail mean/median divergence
# ==============================================================================

cat("\n--- Mean/Median Divergence Check for SAP Bias ---\n")

divergence_table <- results %>%
  group_by(x_type, error_type, outlier_method) %>%
  summarise(
    mean_bias = mean(bias_mis_sap, na.rm = TRUE),
    median_bias = median(bias_mis_sap, na.rm = TRUE),
    ratio = mean_bias / pmax(median_bias, 1e-15),
    n_extreme = sum(
      bias_mis_sap >
        10 * median(bias_mis_sap, na.rm = TRUE),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(ratio > 10 | n_extreme > 5L) %>%
  arrange(desc(ratio))

if (nrow(divergence_table) > 0L) {
  print(divergence_table, n = Inf)
} else {
  cat("No severe mean/median divergence detected.\n")
}


# ==============================================================================
# 14. Save publication tables
# ==============================================================================

safe_save_rds(
  bias_distributional,
  "../output/04sap_bias_distributional.rds"
)

summary_tables <- list(
  health = health_table,
  selected_k = k_table,
  peel_v2_stops = peel_v2_stops,
  sap_stops = sap_stops,
  sap_diagnostics = sap_diagnostics,
  overlap = overlap_table,
  coverage_by_cell = coverage_table,
  coverage_summary = coverage_summary,
  bias_distributional = bias_distributional,
  bias_summary = bias_summary,
  rmse = rmse_table,
  sap_vs_peel = sap_vs_peel,
  runtime = runtime_table,
  sap_bad_leverage = sap_bad_leverage,
  divergence = divergence_table
)

safe_save_rds(
  summary_tables,
  "../output/04sap_summary_tables.rds"
)

cat(
  paste0(
    "\nScript 04 completed successfully.\n",
    "Results: ../output/04sap_robust_comparison_results.rds\n",
    "Tables:  ../output/04sap_summary_tables.rds\n",
    "Bias:    ../output/04sap_bias_distributional.rds\n"
  )
)
