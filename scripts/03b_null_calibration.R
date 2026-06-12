# ==============================================================================
# File: scripts/03b_null_calibration.R
# Purpose: Null calibration test for the EVT test of set influence.
#          Under no contamination, p-values should be Uniform(0,1) and
#          rejection at alpha=0.05 should be ~5%. Tests both constrained
#          and robust GEV fitting strategies across architectures.
#
# Design: For each (architecture, N, k) cell:
#   1. Generate clean data (no injection)
#   2. Pick a random set of size k
#   3. Run EVT (constrained + robust) on that random set
#   4. Record p-values
#   Repeat 500 times per cell.
#
# Outputs:
#   ../output/03b_null_calibration.rds
# ==============================================================================

# 1. Dependencies
library(dplyr)
library(purrr)
library(evd)
library(future)
library(furrr)

source("../R/helpers_local.R")
source("../R/03_scaling_dgp.R")
source("../R/utils_checkpoint.R")
source("../R/exact_dfb_bmx.R")
source("../R/evt_iter.R")
source("../R/evt_iter_dm.R")
source("../R/metrics_influence.R")

# Parallel setup
slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK")
if (slurm_cpus != "") {
  num_workers <- as.numeric(slurm_cpus)
} else {
  num_workers <- max(1, parallel::detectCores() - 2)
}
cat(sprintf("Using %d workers.\n", num_workers))
plan(multisession, workers = num_workers)

# ------------------------------------------------------------------------------
# 2. Configuration
# ------------------------------------------------------------------------------
sim_params <- list(
  n_iters = 500,
  seed    = 20260611
)

set.seed(sim_params$seed)

MIN_OBS_PER_BLOCK <- 30

resolve_block_count <- function(N, k) {
  B_raw <- floor(sqrt(N))
  B_max <- floor((N - k) / MIN_OBS_PER_BLOCK)
  B_max <- max(B_max, 3L)
  min(B_raw, B_max)
}

# Reduced grid — enough to characterize calibration without excess runtime
null_grid <- expand.grid(
  N            = c(500, 2000, 5000),
  k            = c(1, 5, 15),
  architecture = c("simple", "complex", "interaction",
                   "triple_interaction", "nonlinear_nuisance",
                   "sparse_binary_interaction",
                   "polynomial_interaction",
                   "plm_confounded", "plm_nonlinear"),
  stringsAsFactors = FALSE
)

# Filter infeasible cells
null_grid$B <- mapply(resolve_block_count, null_grid$N, null_grid$k)
null_grid <- null_grid[null_grid$B >= 3, ]
rownames(null_grid) <- NULL

n_cells   <- nrow(null_grid)
n_iters   <- sim_params$n_iters
total_its <- n_cells * n_iters

dir.create("../output/temp_03b", recursive = TRUE, showWarnings = FALSE)

cat(sprintf(paste0(
  "Starting 03b Null Calibration.\n",
  "  Grid cells: %d | Iterations per cell: %d | Total: %d\n\n"),
  n_cells, n_iters, total_its))

# ------------------------------------------------------------------------------
# 3. Single Null Iteration Worker
# ------------------------------------------------------------------------------
run_null_iteration <- function(iter_id, N, k, B, architecture, iter_seed) {
  set.seed(iter_seed)
  
  fail_row <- data.frame(
    iter = iter_id,
    p_constrained = NA_real_, conv_constrained = FALSE, shape_constrained = NA_real_,
    p_robust      = NA_real_, conv_robust      = FALSE, shape_robust      = NA_real_,
    error_msg = NA_character_,
    stringsAsFactors = FALSE
  )
  
  tryCatch({
    
    # A. Generate clean data (no injection)
    dgp <- generate_scaling_dgp(N = N, architecture = architecture)
    mod <- lm(dgp$formula, data = dgp$data)
    tpos <- dgp$target_pos
    
    # B. Extract FWL components
    X_full   <- model.matrix(dgp$formula, dgp$data)
    x_target <- X_full[, tpos]
    Z_fwl    <- X_full[, -tpos, drop = FALSE]
    y_vec    <- dgp$data$y
    
    # C. Pick a random set (no outliers — these are ordinary observations)
    random_set <- sample(seq_len(N), k)
    
    # D. EVT with constrained GEV
    evt_con <- tryCatch(
      evt_iter_dm_v2(y = y_vec, x = x_target, Z = Z_fwl,
                     set = random_set, block_count = B,
                     gev_method = "constrained"),
      error = function(e) data.frame(
        shape = NA_real_, scale = NA_real_, loc = NA_real_,
        set_dfb = NA_real_, p_value = NA_real_, converged = FALSE,
        stringsAsFactors = FALSE
      )
    )
    
    # E. EVT with robust GEV
    evt_rob <- tryCatch(
      evt_iter_dm_v2(y = y_vec, x = x_target, Z = Z_fwl,
                     set = random_set, block_count = B,
                     gev_method = "robust"),
      error = function(e) data.frame(
        shape = NA_real_, scale = NA_real_, loc = NA_real_,
        set_dfb = NA_real_, p_value = NA_real_, converged = FALSE,
        stringsAsFactors = FALSE
      )
    )
    
    data.frame(
      iter              = iter_id,
      p_constrained     = evt_con$p_value,
      conv_constrained  = evt_con$converged,
      shape_constrained = evt_con$shape,
      p_robust          = evt_rob$p_value,
      conv_robust       = evt_rob$converged,
      shape_robust      = evt_rob$shape,
      error_msg         = NA_character_,
      stringsAsFactors  = FALSE
    )
    
  }, error = function(e) {
    fail_row$error_msg <- e$message
    fail_row
  })
}

# ------------------------------------------------------------------------------
# 4. Orchestrator Loop
# ------------------------------------------------------------------------------
for (i in seq_len(n_cells)) {
  
  sc <- null_grid[i, ]
  chunk_file <- sprintf("../output/temp_03b/03b_chunk_%03d.rds", i)
  
  if (is_computed(chunk_file)) {
    if (i %% 10 == 0) cat(sprintf("[%03d/%03d] Skipping (cached)\n",
                                  i, n_cells))
    next
  }
  
  cat(sprintf("[%03d/%03d] N=%d k=%d arch=%s B=%d ... ",
              i, n_cells, sc$N, sc$k, sc$architecture, sc$B))
  
  cell_results <- furrr::future_map_dfr(
    seq_len(n_iters), function(iter_id) {
      
      iter_seed <- sim_params$seed + (i - 1) * n_iters + iter_id
      
      run_null_iteration(
        iter_id      = iter_id,
        N            = sc$N,
        k            = sc$k,
        B            = sc$B,
        architecture = sc$architecture,
        iter_seed    = iter_seed
      )
      
    }, .options = furrr_options(seed = TRUE)
  )
  
  cell_results$N            <- sc$N
  cell_results$k            <- sc$k
  cell_results$architecture <- sc$architecture
  cell_results$B            <- sc$B
  
  safe_save_rds(cell_results, chunk_file)
  cat("Done.\n")
}

# ------------------------------------------------------------------------------
# 5. Compile
# ------------------------------------------------------------------------------
cat("\nAssembling null calibration dataset...\n")
null_data <- compile_checkpoints(
  temp_dir          = "../output/temp_03b",
  pattern           = "^03b_chunk_.*\\.rds$",
  final_output_path = "../output/03b_null_calibration.rds",
  clear_temp        = FALSE
)
cat("Compiled successfully.\n")

# ------------------------------------------------------------------------------
# 6. Analysis
# ------------------------------------------------------------------------------
res <- null_data

cat("\n============================\n")
cat("  NULL CALIBRATION RESULTS\n")
cat("============================\n")

cat("\n=== 1. CONVERGENCE ===\n")
conv_table <- res %>%
  group_by(architecture, N) %>%
  summarise(
    conv_con = round(mean(conv_constrained, na.rm = TRUE) * 100, 1),
    conv_rob = round(mean(conv_robust, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )
print(conv_table, n = 50)

cat("\n=== 2. REJECTION RATES AT ALPHA = 0.05 (should be ~5%) ===\n")
rej_table <- res %>%
  group_by(architecture, N, k) %>%
  summarise(
    n_conv_con = sum(conv_constrained, na.rm = TRUE),
    rej_con    = round(mean(p_constrained < 0.05, na.rm = TRUE) * 100, 1),
    n_conv_rob = sum(conv_robust, na.rm = TRUE),
    rej_rob    = round(mean(p_robust < 0.05, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )
cat("Empirical rejection rates (%):\n")
print(rej_table, n = 100)

# Flag problematic cells
overreject <- rej_table %>%
  filter(rej_con > 10 | rej_rob > 10)
if (nrow(overreject) > 0) {
  cat("\nWARNING: Cells with rejection > 10%:\n")
  print(overreject)
}

underreject <- rej_table %>%
  filter((rej_con < 2 & n_conv_con > 50) | (rej_rob < 2 & n_conv_rob > 50))
if (nrow(underreject) > 0) {
  cat("\nWARNING: Cells with rejection < 2% (conservative):\n")
  print(underreject)
}

cat("\n=== 3. KS TEST (p-value uniformity) ===\n")
ks_results <- res %>%
  group_by(architecture, N, k) %>%
  summarise(
    ks_con = {
      pv <- p_constrained[conv_constrained & is.finite(p_constrained)]
      if (length(pv) >= 20) round(ks.test(pv, "punif")$p.value, 3)
      else NA_real_
    },
    ks_rob = {
      pv <- p_robust[conv_robust & is.finite(p_robust)]
      if (length(pv) >= 20) round(ks.test(pv, "punif")$p.value, 3)
      else NA_real_
    },
    .groups = "drop"
  )
cat("KS test p-values (>0.05 = consistent with uniform):\n")
print(ks_results, n = 100)

ks_fail <- ks_results %>%
  filter((!is.na(ks_con) & ks_con < 0.01) | (!is.na(ks_rob) & ks_rob < 0.01))
if (nrow(ks_fail) > 0) {
  cat("\nWARNING: Cells where KS rejects uniformity at 1%:\n")
  print(ks_fail)
} else {
  cat("\nNo cells with severe calibration failure (KS < 0.01).\n")
}

cat("\n=== 4. SHAPE PARAMETERS UNDER NULL ===\n")
shape_null <- res %>%
  filter(conv_constrained == TRUE) %>%
  group_by(architecture) %>%
  summarise(
    shape_con_mean = round(mean(shape_constrained, na.rm = TRUE), 3),
    shape_con_med  = round(median(shape_constrained, na.rm = TRUE), 3),
    shape_rob_mean = round(mean(shape_robust, na.rm = TRUE), 3),
    shape_rob_med  = round(median(shape_robust, na.rm = TRUE), 3),
    .groups = "drop"
  )
cat("Mean/median shape (xi) under null:\n")
print(shape_null)

cat("\n=== 5. OVERALL SUMMARY ===\n")
overall_con <- mean(res$p_constrained < 0.05, na.rm = TRUE) * 100
overall_rob <- mean(res$p_robust < 0.05, na.rm = TRUE) * 100
cat(sprintf("Overall rejection rate: constrained = %.1f%%  robust = %.1f%%\n",
            overall_con, overall_rob))
cat(sprintf("  (Nominal = 5.0%%)\n"))

cat("\nScript 03b complete.\n")