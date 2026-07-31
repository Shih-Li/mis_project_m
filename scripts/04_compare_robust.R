# ==============================================================================
# File: /scripts/04_compare_robust.R
# Purpose: Orchestrate the robust-estimator comparison using exact Dinkelbach
#          MIS detection and Selection-Adjusted Permutation (SAP) calibration.
#
# Estimators (10):
#   full, cd, lev, dfb, mis_alpha, mis_oracle, mis_peel, mis_sap, mm, lts
#
# Outputs:
#   ../output/04_sap_sek_robust_comparison_results.rds
#   ../output/04_sap_sek_summary_tables.rds
#   ../output/04_sap_sek_bias_distributional.rds
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
source("../R/mis_sap.R")
source("../R/mis_sek.R")
source("../R/mis_sek_adapter.R")
source("../R/sim_robust_engine_v2.R")
source("../R/utils_checkpoint.R")

# Guard against accidentally using an old engine left in memory.
required_engine_args <- c(
  "sap_alpha",
  "sap_B_perm",
  "sap_k_grid",
  "sap_max_fraction",
  "sek_alpha",
  "sek_B_cal",
  "sek_k_grid",
  "sek_max_fraction",
  "sek_eta_n",
  "sek_minimum_denominator_fraction"
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
  "iterative_peel_v2",
  "mis_sap",
  "mis_sek",
  "run_mis_sek_from_data",
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
  magnitude = 10,
  seed = 20260503L,
  
  sap_alpha = 0.05,
  sap_B_perm = 199L,
  sap_k_grid = c(1L, 2L, 5L, 10L, 20L, 50L, 100L),
  sap_alpha = 0.05,
  sap_B_perm = 199L,
  sap_k_grid = c(
    1L, 2L, 5L, 10L, 20L, 50L, 100L
  ),
  sap_max_fraction = 0.05,
  
  sek_alpha = 0.05,
  sek_B_cal = 199L,
  sek_k_grid = c(
    1L, 2L, 5L, 10L, 20L, 50L, 100L
  ),
  sek_max_fraction = 0.05,
  sek_eta_n = 0,
  sek_minimum_denominator_fraction = 0.05
)

n_obs_grid <- c(500L, 1000L, 2500L, 5000L)

contam_prop_grid <- c(
  0.005,
  0.010,
  0.025,
  0.050
)

# Grid of N * k
nk_grid <- expand.grid(
  n_obs = n_obs_grid,
  contam_prop = contam_prop_grid,
  stringsAsFactors = FALSE
) %>%
  mutate(
    set_size = pmax(
      floor(n_obs * contam_prop),
      2L
    )
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

contaminated_grid <- merge(
  nk_grid,
  param_grid %>%
    filter(outlier_method != "none"),
  by = NULL
)

clean_grid <- merge(
  data.frame(
    n_obs = n_obs_grid,
    stringsAsFactors = FALSE
  ),
  param_grid %>%
    filter(outlier_method == "none"),
  by = NULL
) %>%
  mutate(
    contam_prop = 0,
    set_size = 0L
  )

design_grid <- bind_rows(
  clean_grid,
  contaminated_grid
) %>%
  arrange(
    n_obs,
    contam_prop,
    x_type,
    error_type,
    outlier_method
  ) %>%
  mutate(
    design_id = row_number()
  )

set.seed(sim_params$seed)

n_designs <- nrow(design_grid)
n_expected <- n_designs * sim_params$n_iters
num_workers <- min(
  max(1L, future::availableCores() - 2L),
  sim_params$n_iters
)

cat(sprintf("Local environment: using %d worker(s).\n", num_workers))
future::plan(future::multisession, workers = num_workers)

cat(sprintf(
  paste0(
    "\nStarting Script 04: Robust MIS-SAP Comparison\n",
    "  Design cells:        %d\n",
    "  Iterations/design:   %d\n",
    "  Expected draws:      %d\n",
    "  Sample sizes:        %s\n",
    "  Contamination grid:  %s\n",
    "  SAP permutations:    %d\n",
    "  SAP k-grid:          %s\n\n"
  ),
  n_designs,
  sim_params$n_iters,
  n_expected,
  paste(n_obs_grid, collapse = ", "),
  paste(
    paste0(100 * contam_prop_grid, "%"),
    collapse = ", "
  ),
  sim_params$sap_B_perm,
  paste(sim_params$sap_k_grid, collapse = ", ")
))


# ==============================================================================
# 3. Checkpointed parallel simulation
# ==============================================================================

for (i in seq_len(n_designs)) {
  d_current <- design_grid[i, , drop = FALSE]
  chunk_file <- sprintf(
    paste0(
      "../output/temp_04_sap_sek/",
      "04_sap_sek_chunk_n%d_cp%04d_k%d_s%04d.rds"
    ),
    d_current$n_obs,
    round(10000 * d_current$contam_prop),
    d_current$set_size,
    i
  )
  
  if (is_computed(chunk_file)) {
    cat(sprintf(
      "[%03d/%03d] Cached: x=%s | error=%s | outlier=%s\n",
      i,
      n_designs,
      d_current$x_type,
      d_current$error_type,
      d_current$outlier_method
    ))
    next
  }
  
  cat(sprintf(
    "[%03d/%03d] Running: x=%s | error=%s | outlier=%s ... ",
    i,
    n_designs,
    d_current$x_type,
    d_current$error_type,
    d_current$outlier_method
  ))
  
  scenario_results <- furrr::future_map_dfr(
    seq_len(sim_params$n_iters),
    function(iter_id) {
      tryCatch(
        run_robust_comparison_iter_v2(
          iter = iter_id,
          n = d_current$n_obs,
          p = 1L,
          x_type = d_current$x_type,
          error_type = d_current$error_type,
          outlier_method = d_current$outlier_method,
          k = d_current$set_size,
          magnitude = sim_params$magnitude,
          sap_alpha = sim_params$sap_alpha,
          sap_B_perm = sim_params$sap_B_perm,
          sap_k_grid =
            sim_params$sap_k_grid,
          
          sap_max_fraction =
            sim_params$sap_max_fraction,
          
          sek_alpha =
            sim_params$sek_alpha,
          
          sek_B_cal =
            sim_params$sek_B_cal,
          
          sek_k_grid =
            sim_params$sek_k_grid,
          
          sek_max_fraction =
            sim_params$sek_max_fraction,
          
          sek_eta_n =
            sim_params$sek_eta_n,
          
          sek_minimum_denominator_fraction =
            sim_params$
            sek_minimum_denominator_fraction
        ),
        error = function(e) {
          warning(sprintf(
            "Iter %d failed for x=%s|error=%s|outlier=%s: %s",
            iter_id,
            d_current$x_type,
            d_current$error_type,
            d_current$outlier_method,
            conditionMessage(e)
          ))
          NULL
        }
      )
    },
    # Scenario-specific seeds remain stable when completed chunks are skipped.
    .options = furrr::furrr_options(seed = sim_params$seed + i)
  )
  
  if (nrow(scenario_results) > 0L) {
    scenario_results <- scenario_results %>%
      mutate(
        n_obs = as.integer(
          d_current$n_obs
        ),
        
        design_k = as.integer(
          d_current$set_size
        ),
        
        # Nominal contamination proportion from the requested grid.
        contam_prop = as.numeric(
          d_current$contam_prop
        ),
        
        # Actual injected proportion after integer rounding of k.
        realized_contam_prop = ifelse(
          outlier_method == "none",
          0,
          set_size / n_obs
        ),
        
        design_id = as.integer(
          d_current$design_id
        )
      )
  }
  
  safe_save_rds(
    scenario_results,
    chunk_file
  )
  
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
  temp_dir = "../output/temp_04_sap_sek",
  pattern = "^04_sap_sek_chunk_.*\\.rds$",
  final_output_path = "../output/04_sap_sek_robust_comparison_results.rds",
  clear_temp = FALSE
)

results <- readRDS("../output/04_sap_sek_robust_comparison_results.rds")

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
  count(
    design_id,
    n_obs,
    design_k,
    contam_prop,
    x_type,
    error_type,
    outlier_method,
    iter,
    name = "n"
  ) %>%
  filter(
    n > 1L
  )

if (nrow(duplicate_rows) > 0L) {
  warning("Duplicate scenario-iteration rows detected.")
  print(duplicate_rows, n = Inf)
}


# ==============================================================================
# 5. Safe summary helpers and health checks
# ==============================================================================

safe_numeric_values <- function(x) {
  x[is.finite(x)]
}


safe_mean <- function(x) {
  x <- safe_numeric_values(x)
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  mean(x)
}


safe_sd <- function(x) {
  x <- safe_numeric_values(x)
  
  if (length(x) < 2L) {
    return(NA_real_)
  }
  
  stats::sd(x)
}


safe_median <- function(x) {
  x <- safe_numeric_values(x)
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  stats::median(x)
}


safe_quantile <- function(x, probability) {
  x <- safe_numeric_values(x)
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  unname(
    stats::quantile(
      x,
      probs = probability,
      na.rm = TRUE,
      names = FALSE
    )
  )
}


safe_rate <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  mean(as.logical(x))
}


safe_rmse <- function(x) {
  x <- safe_numeric_values(x)
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  sqrt(mean(x^2))
}


coverage_label <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0L) {
    return(NA_character_)
  }
  
  p_hat <- mean(x)
  mc_se <- sqrt(
    p_hat *
      (1 - p_hat) /
      length(x)
  )
  
  sprintf(
    "%.1f (%.1f), n=%d",
    100 * p_hat,
    100 * mc_se,
    length(x)
  )
}


cat("\n--- Method Health Checks ---\n")

health_table <- results %>%
  summarise(
    n_rows = n(),
    
    mis_sap_error_rate = safe_rate(
      !mis_sap_ok
    ),
    
    peel_v2_error_rate = safe_rate(
      peel_v2_stop == "error"
    ),
    
    mm_convergence_rate = safe_rate(
      mm_converged
    ),
    
    mis_sap_missing_p_rate = safe_rate(
      !is.finite(mis_sap_global_p)
    ),
    
    mis_sap_missing_coef_rate = safe_rate(
      !is.finite(coef_mis_sap)
    ),
    
    mis_sek_error_rate = safe_rate(
      !mis_sek_ok
    ),
    
    mis_sek_stable_rate = safe_rate(
      mis_sek_state ==
        "stable_effective_k"
    ),
    
    mis_sek_point_rate = safe_rate(
      mis_sek_point_available
    ),
    
    mis_sek_set_rate = safe_rate(
      mis_sek_selection_type == "set"
    ),
    
    mis_sek_invalid_point_coef_rate = {
      point_rows <- which(
        mis_sek_point_available %in% TRUE
      )
      
      if (length(point_rows) == 0L) {
        NA_real_
      } else {
        mean(
          !is.finite(
            coef_mis_sek[point_rows]
          )
        )
      }
    },
    
    mis_sek_unexpected_coef_rate = {
      nonpoint_rows <- which(
        !(mis_sek_point_available %in% TRUE)
      )
      
      if (length(nonpoint_rows) == 0L) {
        NA_real_
      } else {
        mean(
          is.finite(
            coef_mis_sek[nonpoint_rows]
          )
        )
      }
    }
  )

print(health_table)


# ==============================================================================
# 6. Selected-set sizes and process diagnostics
# ==============================================================================

cat("\n--- Selected or Flagged k ---\n")

k_table <- results %>%
  group_by(
    n_obs,
    contam_prop,
    set_size,
    outlier_method
  ) %>%
  summarise(
    true_k = first(set_size),
    
    cd_mean = safe_mean(k_cd),
    lev_mean = safe_mean(k_lev),
    dfb_mean = safe_mean(k_dfb),
    
    alpha_mean = safe_mean(k_alpha),
    oracle_mean = safe_mean(k_oracle),
    
    peel_v2_mean = safe_mean(k_peel_v2),
    peel_v2_med = safe_median(k_peel_v2),
    
    sap_sensitivity_mean = safe_mean(
      k_mis_sap_sensitivity
    ),
    sap_sensitivity_med = safe_median(
      k_mis_sap_sensitivity
    ),
    sap_sensitivity_q25 = safe_quantile(
      k_mis_sap_sensitivity,
      0.25
    ),
    sap_sensitivity_q75 = safe_quantile(
      k_mis_sap_sensitivity,
      0.75
    ),
    
    sap_cleaning_mean = safe_mean(
      k_mis_sap_cleaning
    ),
    sap_cleaning_med = safe_median(
      k_mis_sap_cleaning
    ),
    
    sek_point_k_mean = safe_mean(
      k_mis_sek
    ),
    sek_point_k_med = safe_median(
      k_mis_sek
    ),
    sek_point_rate = safe_rate(
      mis_sek_point_available
    ),
    
    .groups = "drop"
  )

print(k_table, n = Inf)


cat("\n--- Peel-v2 Stop Reasons ---\n")

peel_v2_stops <- results %>%
  count(
    n_obs,
    contam_prop,
    set_size,
    outlier_method,
    peel_v2_stop,
    name = "n"
  ) %>%
  group_by(
    n_obs,
    contam_prop,
    set_size,
    outlier_method
  ) %>%
  mutate(
    pct = round(
      100 * n / sum(n),
      1
    )
  ) %>%
  ungroup() %>%
  arrange(
    n_obs,
    contam_prop,
    outlier_method,
    desc(n)
  )

print(peel_v2_stops, n = Inf)


cat("\n--- MIS-SAP States and Stop Reasons ---\n")

sap_stops <- results %>%
  count(
    n_obs,
    contam_prop,
    set_size,
    outlier_method,
    mis_sap_state,
    mis_sap_stop_reason,
    name = "n"
  ) %>%
  group_by(
    n_obs,
    contam_prop,
    set_size,
    outlier_method
  ) %>%
  mutate(
    pct = round(
      100 * n / sum(n),
      1
    )
  ) %>%
  ungroup() %>%
  arrange(
    n_obs,
    contam_prop,
    outlier_method,
    desc(n)
  )

print(sap_stops, n = Inf)


cat("\n--- MIS-SAP Process Summary ---\n")

sap_diagnostics <- results %>%
  group_by(
    n_obs,
    contam_prop,
    set_size,
    outlier_method
  ) %>%
  summarise(
    detection_rate = safe_rate(
      mis_sap_reject
    ),
    
    cleaning_permitted_rate = safe_rate(
      mis_sap_cleaning_permitted
    ),
    
    exact_k_given_rejection = {
      rejected <- (
        mis_sap_reject %in% TRUE
      )
      
      safe_rate(
        k_mis_sap_sensitivity[rejected] ==
          set_size[rejected]
      )
    },
    
    sensitivity_k_mean = safe_mean(
      k_mis_sap_sensitivity
    ),
    
    sensitivity_k_med = safe_median(
      k_mis_sap_sensitivity
    ),
    
    cleaning_k_mean = safe_mean(
      k_mis_sap_cleaning
    ),
    
    cleaning_k_med = safe_median(
      k_mis_sap_cleaning
    ),
    
    global_p_mean = safe_mean(
      mis_sap_global_p
    ),
    
    global_p_med = safe_median(
      mis_sap_global_p
    ),
    
    peak_excess_mean = safe_mean(
      mis_sap_peak_excess
    ),
    
    peak_excess_med = safe_median(
      mis_sap_peak_excess
    ),
    
    .groups = "drop"
  )

print(sap_diagnostics, n = Inf)


cat("\n--- MIS-sek Certification Summary ---\n")

sek_diagnostics <- results %>%
  group_by(
    n_obs,
    contam_prop,
    set_size,
    outlier_method
  ) %>%
  summarise(
    global_rejection_rate = safe_rate(
      mis_sek_global_reject
    ),
    
    stable_rate = safe_rate(
      mis_sek_state ==
        "stable_effective_k"
    ),
    
    point_identification_rate = safe_rate(
      mis_sek_point_available
    ),
    
    set_identification_rate = safe_rate(
      mis_sek_selection_type == "set"
    ),
    
    support_complete_rate = safe_rate(
      mis_sek_support_complete
    ),
    
    direction_stable_rate = safe_rate(
      mis_sek_direction_stable
    ),
    
    denominator_safe_rate = safe_rate(
      mis_sek_denominator_safe
    ),
    
    formal_guarantee_eligible_rate = safe_rate(
      mis_sek_formal_guarantee_eligible
    ),
    
    point_k_mean = safe_mean(
      mis_sek_selected_k
    ),
    
    point_k_med = safe_median(
      mis_sek_selected_k
    ),
    
    .groups = "drop"
  )

print(sek_diagnostics, n = Inf)


cat("\n--- MIS-sek States ---\n")

sek_states <- results %>%
  count(
    n_obs,
    contam_prop,
    set_size,
    outlier_method,
    mis_sek_state,
    mis_sek_selection_type,
    name = "n"
  ) %>%
  group_by(
    n_obs,
    contam_prop,
    set_size,
    outlier_method
  ) %>%
  mutate(
    pct = round(
      100 * n / sum(n),
      1
    )
  ) %>%
  ungroup() %>%
  arrange(
    n_obs,
    contam_prop,
    outlier_method,
    desc(n)
  )

print(sek_states, n = Inf)


# ==============================================================================
# 7. Detection overlap
# ==============================================================================

cat("\n--- Detection Overlap: Mean (SD) ---\n")

overlap_table <- results %>%
  filter(
    outlier_method != "none"
  ) %>%
  group_by(
    n_obs,
    contam_prop,
    set_size,
    outlier_method
  ) %>%
  summarise(
    cd_mean = safe_mean(overlap_cd),
    cd_sd = safe_sd(overlap_cd),
    
    lev_mean = safe_mean(overlap_lev),
    lev_sd = safe_sd(overlap_lev),
    
    dfb_mean = safe_mean(overlap_dfb),
    dfb_sd = safe_sd(overlap_dfb),
    
    alpha_mean = safe_mean(
      overlap_mis_alpha
    ),
    alpha_sd = safe_sd(
      overlap_mis_alpha
    ),
    
    oracle_mean = safe_mean(
      overlap_mis_oracle
    ),
    oracle_sd = safe_sd(
      overlap_mis_oracle
    ),
    
    peel_v2_mean = safe_mean(
      overlap_peel_v2
    ),
    peel_v2_sd = safe_sd(
      overlap_peel_v2
    ),
    
    sap_sensitivity_mean = safe_mean(
      overlap_mis_sap_sensitivity
    ),
    sap_sensitivity_sd = safe_sd(
      overlap_mis_sap_sensitivity
    ),
    
    sap_cleaning_mean = safe_mean(
      overlap_mis_sap_cleaning
    ),
    sap_cleaning_sd = safe_sd(
      overlap_mis_sap_cleaning
    ),
    
    sek_point_mean = safe_mean(
      overlap_mis_sek
    ),
    sek_point_sd = safe_sd(
      overlap_mis_sek
    ),
    
    .groups = "drop"
  )

print(overlap_table, n = Inf)


# ==============================================================================
# 8. Naive interval coverage
# ==============================================================================

cat(
  paste0(
    "\n--- Naive Post-deletion/Estimator Interval Coverage ",
    "by Simulation Cell (%) ---\n"
  )
)

coverage_table <- results %>%
  group_by(
    n_obs,
    contam_prop,
    set_size,
    x_type,
    error_type,
    outlier_method
  ) %>%
  summarise(
    n_iter = n(),
    
    across(
      starts_with("cov_"),
      ~ {
        valid <- .x[!is.na(.x)]
        
        if (length(valid) == 0L) {
          NA_real_
        } else {
          mean(valid) * 100
        }
      }
    ),
    
    .groups = "drop"
  )

print(coverage_table, n = Inf)


cat(
  paste0(
    "\n--- Coverage by N, k, and Contamination: ",
    "Mean% (Monte Carlo SE%), n ---\n"
  )
)

coverage_summary <- results %>%
  group_by(
    n_obs,
    contam_prop,
    set_size,
    outlier_method
  ) %>%
  summarise(
    across(
      starts_with("cov_"),
      coverage_label
    ),
    
    .groups = "drop"
  )

print(coverage_summary, n = Inf)


# ==============================================================================
# 9. Absolute bias and RMSE
# ==============================================================================

cat("\n--- Absolute Bias by Simulation Cell ---\n")

bias_distributional <- results %>%
  group_by(
    n_obs,
    contam_prop,
    set_size,
    x_type,
    error_type,
    outlier_method
  ) %>%
  summarise(
    across(
      starts_with("bias_"),
      list(
        mean = safe_mean,
        sd = safe_sd,
        med = safe_median,
        q25 = ~ safe_quantile(.x, 0.25),
        q75 = ~ safe_quantile(.x, 0.75)
      ),
      .names = "{.col}__{.fn}"
    ),
    
    .groups = "drop"
  )


cat(
  "\n--- Absolute Bias by N, k, and Contamination: Mean [Median] ---\n"
)

bias_summary <- results %>%
  group_by(
    n_obs,
    contam_prop,
    set_size,
    outlier_method
  ) %>%
  summarise(
    across(
      starts_with("bias_"),
      list(
        mean = safe_mean,
        med = safe_median
      ),
      .names = "{.col}__{.fn}"
    ),
    
    .groups = "drop"
  )

print(bias_summary, n = Inf)


cat("\n--- RMSE by Simulation Cell ---\n")

# Squared absolute bias equals squared estimation error, so this avoids
# hard-coding the true coefficient as 1.
rmse_table <- results %>%
  group_by(
    n_obs,
    contam_prop,
    set_size,
    x_type,
    error_type,
    outlier_method
  ) %>%
  summarise(
    across(
      starts_with("bias_"),
      safe_rmse,
      .names = "rmse_{.col}"
    ),
    
    .groups = "drop"
  ) %>%
  rename_with(
    ~ sub(
      "^rmse_bias_",
      "rmse_",
      .x
    ),
    starts_with("rmse_bias_")
  )

print(rmse_table, n = Inf)


# ==============================================================================
# 10. MIS-SAP versus peel-v2
# ==============================================================================

cat("\n--- MIS-SAP versus Peel-v2 ---\n")

sap_vs_peel <- results %>%
  filter(
    outlier_method != "none"
  ) %>%
  group_by(
    n_obs,
    contam_prop,
    set_size,
    outlier_method
  ) %>%
  summarise(
    k_peel_v2_mean = safe_mean(
      k_peel_v2
    ),
    
    k_peel_v2_med = safe_median(
      k_peel_v2
    ),
    
    k_sap_sensitivity_mean = safe_mean(
      k_mis_sap_sensitivity
    ),
    
    k_sap_sensitivity_med = safe_median(
      k_mis_sap_sensitivity
    ),
    
    k_sap_cleaning_mean = safe_mean(
      k_mis_sap_cleaning
    ),
    
    k_sap_cleaning_med = safe_median(
      k_mis_sap_cleaning
    ),
    
    overlap_peel_v2 = safe_mean(
      overlap_peel_v2
    ),
    
    overlap_sap_sensitivity = safe_mean(
      overlap_mis_sap_sensitivity
    ),
    
    overlap_sap_cleaning = safe_mean(
      overlap_mis_sap_cleaning
    ),
    
    bias_peel_v2 = safe_mean(
      bias_mis_peel
    ),
    
    bias_sap = safe_mean(
      bias_mis_sap
    ),
    
    coverage_peel_v2 = 100 * safe_mean(
      cov_mis_peel
    ),
    
    coverage_sap = 100 * safe_mean(
      cov_mis_sap
    ),
    
    cpu_peel_v2_med = safe_median(
      cpu_peel_v2
    ),
    
    cpu_sap_med = safe_median(
      cpu_mis_sap
    ),
    
    .groups = "drop"
  )

print(sap_vs_peel, n = Inf)


# ==============================================================================
# 11. Complete method runtime
# ==============================================================================

cat("\n--- Complete Runtime: Median [Q25, Q75] Seconds ---\n")

runtime_table <- results %>%
  summarise(
    across(
      starts_with("cpu_"),
      list(
        med = safe_median,
        q25 = ~ safe_quantile(.x, 0.25),
        q75 = ~ safe_quantile(.x, 0.75)
      ),
      .names = "{.col}__{.fn}"
    )
  )

print(
  t(runtime_table)
)


# ==============================================================================
# 12. MIS-SAP and MIS-sek under bad leverage
# ==============================================================================

cat(
  "\n--- MIS-SAP and MIS-sek under Bad Leverage by Error Distribution ---\n"
)

sap_bad_leverage <- results %>%
  filter(
    outlier_method == "bad_leverage"
  ) %>%
  group_by(
    n_obs,
    contam_prop,
    set_size,
    error_type
  ) %>%
  summarise(
    sap_detection_rate = safe_rate(
      mis_sap_reject
    ),
    
    sap_cleaning_rate = safe_rate(
      mis_sap_cleaning_permitted
    ),
    
    sap_exact_k_given_rejection = {
      rejected <- (
        mis_sap_reject %in% TRUE
      )
      
      safe_rate(
        k_mis_sap_sensitivity[rejected] ==
          set_size[rejected]
      )
    },
    
    sap_sensitivity_k_mean = safe_mean(
      k_mis_sap_sensitivity
    ),
    
    sap_sensitivity_k_med = safe_median(
      k_mis_sap_sensitivity
    ),
    
    sap_cleaning_k_mean = safe_mean(
      k_mis_sap_cleaning
    ),
    
    sap_cleaning_k_med = safe_median(
      k_mis_sap_cleaning
    ),
    
    sap_sensitivity_overlap_mean = safe_mean(
      overlap_mis_sap_sensitivity
    ),
    
    sap_cleaning_overlap_mean = safe_mean(
      overlap_mis_sap_cleaning
    ),
    
    sap_bias_mean = safe_mean(
      bias_mis_sap
    ),
    
    sap_bias_med = safe_median(
      bias_mis_sap
    ),
    
    coverage_sap = 100 * safe_mean(
      cov_mis_sap
    ),
    
    coverage_peel_v2 = 100 * safe_mean(
      cov_mis_peel
    ),
    
    coverage_mm = 100 * safe_mean(
      cov_mm
    ),
    
    coverage_oracle = 100 * safe_mean(
      cov_mis_oracle
    ),
    
    sap_global_p_med = safe_median(
      mis_sap_global_p
    ),
    
    sap_peak_excess_med = safe_median(
      mis_sap_peak_excess
    ),
    
    sek_global_rejection_rate = safe_rate(
      mis_sek_global_reject
    ),
    
    sek_stable_rate = safe_rate(
      mis_sek_state ==
        "stable_effective_k"
    ),
    
    sek_point_rate = safe_rate(
      mis_sek_point_available
    ),
    
    sek_set_rate = safe_rate(
      mis_sek_selection_type == "set"
    ),
    
    sek_point_k_mean = safe_mean(
      mis_sek_selected_k
    ),
    
    sek_point_overlap_mean = safe_mean(
      overlap_mis_sek
    ),
    
    .groups = "drop"
  )

print(sap_bad_leverage, n = Inf)


# ==============================================================================
# 13. Heavy-tail mean/median divergence
# ==============================================================================

cat("\n--- Mean/Median Divergence Check for MIS-SAP Bias ---\n")

divergence_table <- results %>%
  group_by(
    n_obs,
    contam_prop,
    set_size,
    x_type,
    error_type,
    outlier_method
  ) %>%
  summarise(
    mean_bias = safe_mean(
      bias_mis_sap
    ),
    
    median_bias = safe_median(
      bias_mis_sap
    ),
    
    ratio = if (
      is.finite(mean_bias) &&
      is.finite(median_bias)
    ) {
      mean_bias /
        pmax(
          median_bias,
          1e-15
        )
    } else {
      NA_real_
    },
    
    n_extreme = {
      valid_bias <- bias_mis_sap[
        is.finite(bias_mis_sap)
      ]
      
      if (length(valid_bias) == 0L) {
        0L
      } else {
        reference_median <- stats::median(
          valid_bias
        )
        
        sum(
          valid_bias >
            10 *
            pmax(
              reference_median,
              1e-15
            )
        )
      }
    },
    
    .groups = "drop"
  ) %>%
  filter(
    ratio > 10 |
      n_extreme > 5L
  ) %>%
  arrange(
    desc(ratio)
  )

if (nrow(divergence_table) > 0L) {
  print(
    divergence_table,
    n = Inf
  )
} else {
  cat(
    "No severe mean/median divergence detected.\n"
  )
}


# ==============================================================================
# 14. Save publication tables
# ==============================================================================

safe_save_rds(
  bias_distributional,
  "../output/04_sap_sek_bias_distributional.rds"
)

summary_tables <- list(
  health = health_table,
  
  selected_k = k_table,
  
  peel_v2_stops = peel_v2_stops,
  
  sap_stops = sap_stops,
  sap_diagnostics = sap_diagnostics,
  
  sek_states = sek_states,
  sek_diagnostics = sek_diagnostics,
  
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
  "../output/04_sap_sek_summary_tables.rds"
)

cat(
  paste0(
    "\nScript 04 completed successfully.\n",
    "Results: ../output/04_sap_sek_robust_comparison_results.rds\n",
    "Tables:  ../output/04_sap_sek_summary_tables.rds\n",
    "Bias:    ../output/04_sap_sek_bias_distributional.rds\n"
  )
)