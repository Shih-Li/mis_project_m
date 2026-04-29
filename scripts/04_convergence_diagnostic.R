# ==============================================================================
# File:    script/03_dryrun_test.R
# Purpose: Verify the B-cap fix before committing to a full rerun.
#          Part 1: Pure arithmetic — no simulation, no sourcing, instant.
#          Part 2: Run 10 iterations on the 10 worst-case cells from diagnostic.
#          Total runtime: ~5 min.
#
# Usage:   Run from the script/ directory (same working directory as 03_alg_comp.R).
#          Requires the same ../R/ source files and packages as 03_alg_comp.R.
# ==============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("  03 DRY RUN TEST\n")
cat(strrep("=", 70), "\n\n")

# ==============================================================================
# PART 1: resolve_block_count arithmetic check (instant, no dependencies)
# ==============================================================================

cat("--- PART 1: B-cap arithmetic check ---\n\n")

MIN_OBS_PER_BLOCK <- 30

resolve_block_count <- function(B_type, N, k) {
  B_raw <- if (B_type == "sqrt") floor(sqrt(N)) else as.numeric(B_type)
  B_max <- floor((N - k) / MIN_OBS_PER_BLOCK)
  B_max <- max(B_max, 3L)
  min(B_raw, B_max)
}

# Rebuild the full scenario grid (same logic as 03_alg_comp.R section 5)
grid_base <- expand.grid(
  N = c(500, 1000, 2000, 5000),
  k = c(1, 3, 5, 10, 15, 20),
  B_type = c("20", "50", "100", "sqrt"),
  architecture = c("simple","complex","interaction","triple_interaction",
                   "nonlinear_nuisance","sparse_binary_interaction",
                   "polynomial_interaction"),
  stringsAsFactors = FALSE
)

grid_high_k <- expand.grid(
  N = c(500, 1000, 2000, 5000),
  k = c(1, 3, 5, 10, 15, 20, 30, 50),
  B_type = c("20", "50", "100", "sqrt"),
  architecture = "high_k_interaction",
  stringsAsFactors = FALSE
)

grid_collinear <- expand.grid(
  N = c(500, 1000, 2000, 5000),
  k = c(1, 3, 5, 10, 15, 20),
  B_type = c("20", "50", "100", "sqrt"),
  architecture = "collinear_interaction",
  rho = c(0.5, 0.7, 0.85, 0.95),
  stringsAsFactors = FALSE
)

grid_base$rho    <- NA_real_
grid_high_k$rho  <- NA_real_
scenario_grid <- rbind(grid_base, grid_high_k, grid_collinear)
scenario_grid <- scenario_grid[scenario_grid$k / scenario_grid$N <= 0.05, ]
rownames(scenario_grid) <- NULL

# Compute B_raw and B_actual for every scenario
scenario_grid$B_raw <- ifelse(
  scenario_grid$B_type == "sqrt",
  floor(sqrt(scenario_grid$N)),
  as.numeric(scenario_grid$B_type)
)
scenario_grid$B_actual <- mapply(
  resolve_block_count,
  scenario_grid$B_type, scenario_grid$N, scenario_grid$k
)
scenario_grid$B_capped <- scenario_grid$B_actual < scenario_grid$B_raw
scenario_grid$obs_per_block <- floor(
  (scenario_grid$N - scenario_grid$k) / scenario_grid$B_actual
)

# ── Check 1: No scenario has fewer than 3 blocks ──
bad_blocks <- scenario_grid[scenario_grid$B_actual < 3, ]
if (nrow(bad_blocks) > 0) {
  cat("[FAIL] Scenarios with B_actual < 3:\n")
  print(bad_blocks[, c("N","k","B_type","B_raw","B_actual")], row.names = FALSE)
} else {
  cat("[PASS] All scenarios have B_actual >= 3.\n")
}

# ── Check 2: Every scenario has >= MIN_OBS_PER_BLOCK obs per block ──
# (except where B_raw was already <= 3, in which case capping can't help)
thin <- scenario_grid[scenario_grid$obs_per_block < MIN_OBS_PER_BLOCK, ]
if (nrow(thin) > 0) {
  cat(sprintf("[WARN] %d scenarios still have < %d obs/block after capping:\n",
              nrow(thin), MIN_OBS_PER_BLOCK))
  print(
    unique(thin[, c("N","k","B_type","B_raw","B_actual","obs_per_block")]),
    row.names = FALSE
  )
} else {
  cat(sprintf("[PASS] All scenarios have >= %d obs/block after capping.\n",
              MIN_OBS_PER_BLOCK))
}

# ── Check 3: Summary of capping ──
n_capped <- sum(scenario_grid$B_capped)
cat(sprintf("\nCapping summary: %d / %d scenarios capped (%.1f%%)\n",
            n_capped, nrow(scenario_grid),
            n_capped / nrow(scenario_grid) * 100))

# Show unique capped (N, k, B_type) combos
if (n_capped > 0) {
  capped <- unique(scenario_grid[scenario_grid$B_capped,
                                 c("N","k","B_type","B_raw","B_actual","obs_per_block")])
  capped <- capped[order(capped$N, capped$k, capped$B_type), ]
  cat("\nCapped combinations:\n")
  print(capped, row.names = FALSE)
}

# ── Check 4: The specific worst cells from diagnostic ──
cat("\n--- Spot-check: diagnostic 04 worst cells ---\n")
worst_cases <- data.frame(
  N      = c(500,  500,  500,  500, 1000, 1000),
  k      = c(5,    10,   15,   20,  10,   15),
  B_type = c("100","100","100","100","100","100"),
  label  = c("was 0%","was 0%","was 0%","was 0%","was 0%","was 0%"),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(worst_cases))) {
  w <- worst_cases[i, ]
  B <- resolve_block_count(w$B_type, w$N, w$k)
  B_raw <- as.numeric(w$B_type)
  opb <- floor((w$N - w$k) / B)
  cat(sprintf("  N=%d k=%d B_type=%s: B_raw=%d -> B_actual=%d (%d obs/block) %s\n",
              w$N, w$k, w$B_type, B_raw, B, opb, w$label))
}

cat("\n[PART 1 COMPLETE]\n")

# ==============================================================================
# PART 2: Live test on worst-case cells (10 iterations each)
# ==============================================================================

cat("\n--- PART 2: Live simulation test (10 worst cells x 10 iters) ---\n\n")

library(dplyr)
library(purrr)
library(testingMIS)
library(influence)
library(evd)
library(future)
library(furrr)

source("../R/03_scaling_dgp.R")
source("../R/utils_checkpoint.R")
source("../R/exact_dfb_bmx.R")
source("../R/evt_iter_dm.R")
source("../R/evt_iter.R")

plan(multisession, workers = 14)

# Copy inject_safe_outliers from 03_alg_comp.R
inject_safe_outliers <- function(dgp_res, k, magnitude = 4) {
  df <- dgp_res$data; N <- nrow(df)
  outlier_idx <- sample(1:N, k, replace = FALSE)
  shift_x <- max(mad(df$X, constant = 1.4826), 1e-4) * magnitude
  df$X[outlier_idx] <- df$X[outlier_idx] + shift_x
  if ("M" %in% names(df)) {
    if (length(unique(df$M)) <= 2) {
      df$M[outlier_idx] <- 1L
    } else {
      shift_m <- max(mad(df$M, constant = 1.4826), 1e-4) * magnitude
      df$M[outlier_idx] <- df$M[outlier_idx] + shift_m
    }
  }
  if ("W" %in% names(df)) {
    shift_w <- max(mad(df$W, constant = 1.4826), 1e-4) * magnitude
    df$W[outlier_idx] <- df$W[outlier_idx] + shift_w
  }
  z_cols <- grep("^(Z_matrix\\.|X[0-9]+|Z[0-9]+)$", names(df), value = TRUE)
  if (length(z_cols) > 0) {
    for (zc in head(z_cols, 3)) {
      shift_z <- max(mad(df[[zc]], constant = 1.4826), 1e-4) * magnitude
      df[[zc]][outlier_idx] <- df[[zc]][outlier_idx] + shift_z
    }
  }
  clean_mod  <- lm(dgp_res$formula, data = df[-outlier_idx, ])
  expected_y <- predict(clean_mod, newdata = df[outlier_idx, ])
  scale_y    <- max(mad(df$y, constant = 1.4826), 1e-4)
  df$y[outlier_idx] <- expected_y - (sign(expected_y) * scale_y * magnitude * 2)
  dgp_res$data <- df; dgp_res$true_outliers <- outlier_idx
  return(dgp_res)
}

# Copy run_scaling_iteration from 03_alg_comp.R (verbatim)
run_scaling_iteration <- function(iter_id, N, k, B, architecture,
                                  magnitude, iter_seed, rho = NA) {
  set.seed(iter_seed)
  tryCatch({
    dgp_clean <- generate_scaling_dgp(
      N = N, architecture = architecture,
      rho = if (!is.null(rho) && !is.na(rho)) rho else 0.8
    )
    dgp_poisoned <- inject_safe_outliers(dgp_clean, k = k, magnitude = magnitude)
    mod  <- lm(dgp_poisoned$formula, data = dgp_poisoned$data)
    tpos <- dgp_poisoned$target_pos
    sens_obj <- influence::sens(
      mod, lambda = influence::set_lambda("beta_i", pos = tpos,
                                          sign = sign(coef(mod)[tpos])))
    detected_set <- sens_obj$influence$id[1:k]
    X_full   <- model.matrix(dgp_poisoned$formula, dgp_poisoned$data)
    x_target <- X_full[, tpos]
    Z_fwl    <- X_full[, -tpos, drop = FALSE]
    y_vec    <- dgp_poisoned$data$y
    
    t_greedy_start <- Sys.time()
    res_greedy <- evt_iter(y = y_vec, x = x_target, Z = Z_fwl,
                           set = detected_set, block_count = B)
    t_greedy <- as.numeric(difftime(Sys.time(), t_greedy_start, units = "secs"))
    
    sens_pos <- influence::sens(mod, lambda = influence::set_lambda(
      "beta_i", pos = tpos, sign = 1))
    sens_neg <- influence::sens(mod, lambda = influence::set_lambda(
      "beta_i", pos = tpos, sign = -1))
    detected_pos <- sens_pos$influence$id[1:k]
    detected_neg <- sens_neg$influence$id[1:k]
    overlap_pos <- length(intersect(detected_pos, dgp_poisoned$true_outliers))
    overlap_neg <- length(intersect(detected_neg, dgp_poisoned$true_outliers))
    detected_set_exact <- if (overlap_neg > overlap_pos) detected_neg else detected_pos
    
    t_exact_start <- Sys.time()
    res_exact <- evt_iter_dm(y = y_vec, x = x_target, Z = Z_fwl,
                             set = detected_set_exact, block_count = B)
    t_exact <- as.numeric(difftime(Sys.time(), t_exact_start, units = "secs"))
    
    det_greedy <- length(intersect(detected_set, dgp_poisoned$true_outliers)) / k
    det_exact  <- length(intersect(detected_set_exact, dgp_poisoned$true_outliers)) / k
    
    data.frame(
      iter               = iter_id,
      detection_rate       = det_greedy,
      detection_rate_exact = det_exact,
      p_greedy         = res_greedy$p_value,
      converged_greedy = res_greedy$converged,
      p_exact          = res_exact$p_value,
      converged_exact  = res_exact$converged,
      cpu_greedy       = t_greedy,
      cpu_exact        = t_exact,
      error_msg        = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning(sprintf("Iter %d failed: %s", iter_id, e$message))
    data.frame(
      iter               = iter_id,
      detection_rate       = NA_real_,
      detection_rate_exact = NA_real_,
      p_greedy         = NA_real_,
      converged_greedy = FALSE,
      p_exact          = NA_real_,
      converged_exact  = FALSE,
      cpu_greedy       = NA_real_,
      cpu_exact        = NA_real_,
      error_msg        = e$message,
      stringsAsFactors = FALSE
    )
  })
}

# ── Define 10 worst-case test cells ──
# These all had 0% exact convergence in diagnostic Phase 1.
# Mix of architectures so we're not just testing 'simple'.
test_cells <- data.frame(
  N            = c( 500,  500,  500,  500, 1000, 1000,  500,  500,  500, 1000),
  k            = c(   5,   10,   15,   20,   10,   15,   10,   15,   20,   20),
  B_type       = c("100","100","100","100","100","100", "50", "50", "50", "50"),
  architecture = c("simple","simple","simple","simple","simple","simple",
                   "simple","simple","simple","simple"),
  rho          = NA_real_,
  stringsAsFactors = FALSE
)

N_TEST_ITERS <- 10
SEED_BASE    <- 99990000L

cat(sprintf("Running %d test cells x %d iterations = %d runs...\n\n",
            nrow(test_cells), N_TEST_ITERS,
            nrow(test_cells) * N_TEST_ITERS))

all_results <- list()
t0 <- proc.time()

for (tc in seq_len(nrow(test_cells))) {
  cell <- test_cells[tc, ]
  
  # Old B (uncapped)
  B_old <- if (cell$B_type == "sqrt") floor(sqrt(cell$N)) else as.numeric(cell$B_type)
  # New B (capped)
  B_new <- resolve_block_count(cell$B_type, cell$N, cell$k)
  was_capped <- B_new < B_old
  
  cat(sprintf("  Cell %2d: N=%d k=%d B_type=%s  B_old=%d  B_new=%d %s\n",
              tc, cell$N, cell$k, cell$B_type, B_old, B_new,
              if (was_capped) "[CAPPED]" else ""))
  
  # Run with OLD B (uncapped) — the broken behaviour
  res_old <- furrr::future_map_dfr(
    seq_len(N_TEST_ITERS), function(j) {
      run_scaling_iteration(
        iter_id = j, N = cell$N, k = cell$k, B = B_old,
        architecture = cell$architecture, magnitude = 5,
        iter_seed = SEED_BASE + (tc - 1) * 100L + j, rho = cell$rho
      )
    }, .options = furrr_options(seed = TRUE)
  )
  
  # Run with NEW B (capped) — the fix
  res_new <- furrr::future_map_dfr(
    seq_len(N_TEST_ITERS), function(j) {
      run_scaling_iteration(
        iter_id = j, N = cell$N, k = cell$k, B = B_new,
        architecture = cell$architecture, magnitude = 5,
        iter_seed = SEED_BASE + (tc - 1) * 100L + j, rho = cell$rho
      )
    }, .options = furrr_options(seed = TRUE)
  )
  
  conv_old_exact  <- mean(res_old$converged_exact,  na.rm = TRUE)
  conv_new_exact  <- mean(res_new$converged_exact,  na.rm = TRUE)
  conv_old_greedy <- mean(res_old$converged_greedy, na.rm = TRUE)
  conv_new_greedy <- mean(res_new$converged_greedy, na.rm = TRUE)
  det_new         <- mean(res_new$detection_rate_exact, na.rm = TRUE)
  
  all_results[[tc]] <- data.frame(
    N = cell$N, k = cell$k, B_type = cell$B_type,
    B_old = B_old, B_new = B_new, capped = was_capped,
    conv_exact_old   = round(conv_old_exact,  2),
    conv_exact_new   = round(conv_new_exact,  2),
    conv_greedy_old  = round(conv_old_greedy, 2),
    conv_greedy_new  = round(conv_new_greedy, 2),
    det_exact_new    = round(det_new, 3),
    errors_old       = sum(!is.na(res_old$error_msg)),
    errors_new       = sum(!is.na(res_new$error_msg)),
    stringsAsFactors = FALSE
  )
}

plan(sequential)
elapsed <- (proc.time() - t0)[["elapsed"]]

results_df <- do.call(rbind, all_results)

cat("\n")
cat(strrep("=", 70), "\n")
cat("  DRY RUN RESULTS\n")
cat(strrep("=", 70), "\n\n")

print(results_df, row.names = FALSE)

cat(sprintf("\nTotal time: %.1f min\n", elapsed / 60))

# ── Verdict ──
cat("\n--- VERDICT ---\n")

improved <- results_df$conv_exact_new > results_df$conv_exact_old
still_zero <- results_df$conv_exact_new == 0
any_errors <- results_df$errors_new > 0

if (all(improved | results_df$conv_exact_old == results_df$conv_exact_new)) {
  cat("[PASS] Exact convergence improved (or unchanged) in all cells.\n")
} else {
  cat("[WARN] Some cells got worse:\n")
  print(results_df[!improved & results_df$conv_exact_old != results_df$conv_exact_new, ],
        row.names = FALSE)
}

if (any(still_zero)) {
  cat("[WARN] Some cells still at 0% convergence after fix:\n")
  print(results_df[still_zero, c("N","k","B_type","B_old","B_new","conv_exact_new")],
        row.names = FALSE)
} else {
  cat("[PASS] No cells at 0% convergence after fix.\n")
}

if (any(any_errors)) {
  cat("[WARN] Pipeline errors in new runs:\n")
  print(results_df[any_errors, c("N","k","B_type","errors_new")], row.names = FALSE)
} else {
  cat("[PASS] No pipeline errors.\n")
}

# Check detection rate wasn't harmed
low_det <- results_df$det_exact_new < 0.90
if (any(low_det, na.rm = TRUE)) {
  cat("[WARN] Detection rate < 90% in some cells:\n")
  print(results_df[low_det, c("N","k","B_type","det_exact_new")], row.names = FALSE)
} else {
  cat("[PASS] Detection rate >= 90% in all cells.\n")
}

cat("\nIf all PASS: safe to delete ../output/temp_03/ and run 03_alg_comp.R.\n")
cat("If any WARN: investigate before full rerun.\n")