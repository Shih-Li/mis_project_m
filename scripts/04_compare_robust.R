# ==============================================================================
# File: /scripts/04_compare_robust.R
# Purpose:
#   Robust-estimator comparison for Script 04.
#
# Methods:
#   1. OLS
#   2. Leverage deletion + OLS refit
#   3. Cook's-distance deletion + OLS refit
#   4. DFBETAS deletion + OLS refit
#   5. MIS with oracle k and oracle direction + OLS refit
#   6. MM
#   7. LTS
#
# Final output:
#   ../output/04_robust_comparison_results.rds
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
source("../R/dinkelbach_topk.R")
source("../R/sim_robust_engine.R")
source("../R/utils_checkpoint.R")


required_functions <- c(
  "generate_complex_data",
  "apply_influence_shift",
  "get_leverage",
  "get_cooks_d",
  "get_dfbetas",
  "fit_clean_ols",
  "fit_mm_estimator",
  "fit_lts_estimator",
  "dinkelbach_topk_lm",
  "extract_target_lm",
  "check_coverage",
  "run_robust_comparison_iter",
  "safe_save_rds",
  "is_computed",
  "compile_checkpoints"
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
  n_iters = 100L,
  magnitude = 10,
  seed = 20260503L
)

checkpoint_dir <-
  "../output/temp_04_robust_final_signfix"

final_output_path <-
  "../output/04_robust_comparison_results.rds"

dir.create(
  checkpoint_dir,
  recursive = TRUE,
  showWarnings = FALSE
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

if (nrow(design_grid) != 1248L) {
  stop(
    "Expected 1,248 design cells; found ",
    nrow(design_grid),
    "."
  )
}

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
    "\nStarting Script 04: Robust Comparison\n",
    "  Design cells:        %d\n",
    "  Iterations/design:   %d\n",
    "  Expected draws:      %d\n",
    "  Sample sizes:        %s\n",
    "  Contamination grid:  %s\n",
    "  Workers:             %d\n\n"
  ),
  n_designs,
  sim_params$n_iters,
  n_expected,
  paste(n_obs_grid, collapse = ", "),
  paste(
    paste0(
      100 * contam_prop_grid,
      "%"
    ),
    collapse = ", "
  ),
  num_workers
))


# ==============================================================================
# 3. Checkpointed parallel simulation
# ==============================================================================

for (i in seq_len(n_designs)) {
  d_current <- design_grid[i, , drop = FALSE]
  chunk_file <- file.path(
    checkpoint_dir,
    sprintf(
      "04_robust_chunk_n%d_cp%04d_k%d_s%04d.rds",
      d_current$n_obs,
      round(10000 * d_current$contam_prop),
      d_current$set_size,
      i
    )
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
        run_robust_comparison_iter(
          iter = iter_id,
          n = d_current$n_obs,
          p = 1L,
          x_type = d_current$x_type,
          error_type = d_current$error_type,
          outlier_method = d_current$outlier_method,
          k = d_current$set_size,
          magnitude = sim_params$magnitude
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
  
  if (nrow(scenario_results) != sim_params$n_iters) {
    
    warning(sprintf(
      paste0(
        "Design %d produced only %d/%d successful iterations. ",
        "This cell will NOT be checkpointed and will be rerun ",
        "on the next execution."
      ),
      i,
      nrow(scenario_results),
      sim_params$n_iters
    ))
    
    next
  }
  
  scenario_results <- scenario_results %>%
    mutate(
      n_obs = as.integer(
        d_current$n_obs
      ),
      
      design_k = as.integer(
        d_current$set_size
      ),
      
      contam_prop = as.numeric(
        d_current$contam_prop
      ),
      
      realized_contam_prop = ifelse(
        outlier_method == "none",
        0,
        set_size / n_obs
      ),
      
      design_id = as.integer(
        d_current$design_id
      )
    )
  
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
  temp_dir = checkpoint_dir,
  pattern = "^04_robust_chunk_.*\\.rds$",
  final_output_path = final_output_path,
  clear_temp = FALSE
)

results <- readRDS(final_output_path)

cell_counts <- results %>%
  count(
    design_id,
    n_obs,
    design_k,
    contam_prop,
    x_type,
    error_type,
    outlier_method,
    name = "n_iter"
  )

if (nrow(cell_counts) != 1248L) {
  stop(
    "Expected 1,248 completed design cells; found ",
    nrow(cell_counts),
    "."
  )
}

if (any(cell_counts$n_iter != sim_params$n_iters)) {
  stop(
    "At least one design cell does not contain exactly ",
    sim_params$n_iters,
    " iterations."
  )
}

if (nrow(results) != 124800L) {
  stop(
    "Formal Script 04 run must contain exactly 124,800 rows; found ",
    nrow(results),
    "."
  )
}

cat(sprintf(
  "Observed rows: %d | Expected rows: %d | Missing rows: %d\n",
  nrow(results),
  n_expected,
  n_expected - nrow(results)
))

if (nrow(results) != n_expected) {
  stop(
    "Formal Script 04 dataset is incomplete: expected ",
    n_expected,
    " rows but found ",
    nrow(results),
    "."
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
# 5. Method-specific integrity checks
# ==============================================================================

clean_rows <- results$outlier_method == "none"
contam_rows <- !clean_rows

# Clean cells: oracle MIS must equal OLS exactly.
if (!all(results$k_oracle[clean_rows] == 0L)) {
  stop(
    "Clean-data integrity failure: ",
    "k_oracle is not always zero."
  )
}

if (!all(results$oracle_direction[clean_rows] == 0L)) {
  stop(
    "Clean-data integrity failure: ",
    "oracle_direction is not always zero."
  )
}

if (!identical(
  results$coef_mis_oracle[clean_rows],
  results$coef_full[clean_rows]
)) {
  stop(
    "Clean-data integrity failure: ",
    "MIS oracle coefficient does not exactly equal OLS."
  )
}

# Contaminated cells: oracle k must equal injected k.
if (!all(
  results$k_oracle[contam_rows] ==
  results$set_size[contam_rows]
)) {
  stop(
    "Contaminated-data integrity failure: ",
    "oracle k differs from injected set size."
  )
}

# Oracle direction must be valid.
if (!all(
  results$oracle_direction %in%
  c(-1L, 0L, 1L)
)) {
  stop(
    "Invalid oracle_direction detected."
  )
}

# Oracle direction must agree with the Dinkelbach
# full-minus-deleted sign convention.
expected_oracle_direction <- ifelse(
  results$oracle_dfbeta_delta[contam_rows] >= 0,
  1L,
  -1L
)

if (!all(
  results$oracle_direction[contam_rows] ==
  expected_oracle_direction
)) {
  stop(
    "Oracle-direction integrity failure: ",
    "recorded sign does not agree with ",
    "beta_full - beta_deleted."
  )
}

# Clean rows should have no directional shift.
if (!all(
  results$oracle_dfbeta_delta[clean_rows] == 0
)) {
  stop(
    "Clean-data integrity failure: ",
    "oracle_dfbeta_delta is not zero."
  )
}

# No obsolete Script-04 fields should survive.
obsolete_cols <- grep(
  "mis_alpha|sap|sek|peel",
  names(results),
  value = TRUE,
  ignore.case = TRUE
)

if (length(obsolete_cols) > 0L) {
  stop(
    "Obsolete result columns detected: ",
    paste(
      obsolete_cols,
      collapse = ", "
    )
  )
}

cat(
  "\nScript 04 integrity checks passed.\n"
)