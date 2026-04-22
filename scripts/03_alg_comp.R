# ==============================================================================
# File: script/03_alg_comp.R
# Purpose: Execute dimensional scaling and structural robustness simulations.
#          Compares the Greedy heuristic (testingMIS::dfb_bmx) against
#          Exact Dinkelbach (exact_dfb_bmx) across sample sizes, set sizes,
#          block-count strategies, and five model architectures.
#          Quantifies detection accuracy, EVD convergence, and runtime scaling.
# Inputs:  ../R/03_scaling_dgp.R, ../R/exact_dfb_bmx.R, ../R/evt_iter.R,
#          ../R/evt_iter_dm.R, ../R/utils_checkpoint.R
# Outputs: ../output/temp_03_03/03_chunk_*.rds -> ../output/03_scaling_results_master.rds
# Paper Section: Algorithmic Exactness and Computational Scalability
# ==============================================================================

# 1. Load Dependencies & Source Engines
library(dplyr)
library(purrr)
library(testingMIS)
library(influence)
library(evd)
library(future)
library(furrr)

# Source order matters: exact_dfb_bmx BEFORE evt_iter_dm (handoff rule #4)
source("../R/03_scaling_dgp.R")
source("../R/utils_checkpoint.R")
source("../R/exact_dfb_bmx.R")
source("../R/evt_iter_dm.R")
source("../R/evt_iter.R")



# Set up parallel processing (leave 2 threads for OS stability)
plan(multisession, workers = 14)
cat("Parallel processing initialized with 14 workers.\n\n")

# ------------------------------------------------------------------------------
# 2. Global Configuration
# ------------------------------------------------------------------------------
sim_params <- list(
  n_iters   = 1,
  magnitude = 5,
  seed      = 20260421
)

set.seed(sim_params$seed)

# ------------------------------------------------------------------------------
# 3. The Safe Injection Adapter (Preserves Interaction Geometry)
# ------------------------------------------------------------------------------
inject_safe_outliers <- function(dgp_res, k, magnitude = 4) {
  df <- dgp_res$data
  N  <- nrow(df)
  outlier_idx <- sample(1:N, k, replace = FALSE)
  
  # 3a. Shift main covariates robustly
  shift_x <- max(mad(df$X, constant = 1.4826), 1e-4) * magnitude
  df$X[outlier_idx] <- df$X[outlier_idx] + shift_x
  
  if ("M" %in% names(df)) {
    shift_m <- max(mad(df$M, constant = 1.4826), 1e-4) * magnitude
    df$M[outlier_idx] <- df$M[outlier_idx] + shift_m
  }
  if ("W" %in% names(df)) {
    shift_w <- max(mad(df$W, constant = 1.4826), 1e-4) * magnitude
    df$W[outlier_idx] <- df$W[outlier_idx] + shift_w
  }
  
  # 3b. For complex architecture, shift a few Z columns too so outliers are
  #     actually influential in the high-dimensional covariate space
  z_cols <- grep("^(Z_matrix\\.|X[0-9]+|Z[0-9]+)$", names(df), value = TRUE)
  if (length(z_cols) > 0) {
    for (zc in head(z_cols, 3)) {
      shift_z <- max(mad(df[[zc]], constant = 1.4826), 1e-4) * magnitude
      df[[zc]][outlier_idx] <- df[[zc]][outlier_idx] + shift_z
    }
  }
  
  # 3c. Recalculate Y to forcefully mask the target parameter
  clean_mod  <- lm(dgp_res$formula, data = df[-outlier_idx, ])
  expected_y <- predict(clean_mod, newdata = df[outlier_idx, ])
  scale_y    <- max(mad(df$y, constant = 1.4826), 1e-4)
  df$y[outlier_idx] <- expected_y - (sign(expected_y) * scale_y * magnitude * 2)
  
  dgp_res$data <- df
  dgp_res$true_outliers <- outlier_idx
  return(dgp_res)
}

# ------------------------------------------------------------------------------
# 4. Single-Iteration Worker (Called Inside future_map_dfr)
# ------------------------------------------------------------------------------
run_scaling_iteration <- function(iter_id, N, k, B, architecture,
                                  magnitude, iter_seed) {
  set.seed(iter_seed)
  
  tryCatch({
    
    # A. Generate Data
    dgp_clean <- generate_scaling_dgp(N = N, architecture = architecture)
    
    # B. Inject Adversarial Outliers
    dgp_poisoned <- inject_safe_outliers(dgp_clean, k = k, magnitude = magnitude)
    
    # C. Fit base model on poisoned data
    mod  <- lm(dgp_poisoned$formula, data = dgp_poisoned$data)
    tpos <- dgp_poisoned$target_pos
    
    # D. Detect the influential set via influence::sens()
    sens_obj <- influence::sens(
      mod,
      lambda = influence::set_lambda(
        "beta_i", pos = tpos, sign = sign(coef(mod)[tpos])
      )
    )
    detected_set <- sens_obj$influence$id[1:k]
    
    # E. Extract FWL components: x = target column, Z = everything else
    #    Z must NEVER contain x (handoff rule #1)
    X_full   <- model.matrix(dgp_poisoned$formula, dgp_poisoned$data)
    x_target <- X_full[, tpos]
    Z_fwl    <- X_full[, -tpos, drop = FALSE]
    y_vec    <- dgp_poisoned$data$y
    
    # F. Run Greedy EVT
    t_greedy_start <- Sys.time()
    res_greedy <- evt_iter(
      y = y_vec, x = x_target, Z = Z_fwl,
      set = detected_set, block_count = B
    )
    t_greedy <- as.numeric(difftime(Sys.time(), t_greedy_start, units = "secs"))
    
    # G. Dual-sign detection (Oracle: pick the direction that finds more outliers)
    sens_pos <- influence::sens(mod, lambda = influence::set_lambda("beta_i", pos = tpos, sign = 1))
    sens_neg <- influence::sens(mod, lambda = influence::set_lambda("beta_i", pos = tpos, sign = -1))
    
    detected_pos <- sens_pos$influence$id[1:k]
    detected_neg <- sens_neg$influence$id[1:k]
    
    overlap_pos <- length(intersect(detected_pos, dgp_poisoned$true_outliers))
    overlap_neg <- length(intersect(detected_neg, dgp_poisoned$true_outliers))
    
    if (overlap_neg > overlap_pos) {
      detected_set_exact <- detected_neg
    } else {
      detected_set_exact <- detected_pos
    }
    
    # H. Run Exact Dinkelbach EVT on the best-direction set
    t_exact_start <- Sys.time()
    res_exact <- evt_iter_dm(
      y = y_vec, x = x_target, Z = Z_fwl,
      set = detected_set_exact, block_count = B
    )
    t_exact <- as.numeric(difftime(Sys.time(), t_exact_start, units = "secs"))
    
    # I. Detection metrics (greedy = single-sign, exact = oracle dual-sign)
    detection_rate       <- length(intersect(detected_set, dgp_poisoned$true_outliers)) / k
    detection_rate_exact <- length(intersect(detected_set_exact, dgp_poisoned$true_outliers)) / k
    
    data.frame(
      iter               = iter_id,
      detection_rate       = detection_rate,
      detection_rate_exact = detection_rate_exact,
      p_greedy         = res_greedy$p_value,
      converged_greedy = res_greedy$converged,
      p_exact          = res_exact$p_value,
      converged_exact  = res_exact$converged,
      cpu_greedy   = t_greedy,
      cpu_exact    = t_exact,
      error_msg    = NA_character_,
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
      cpu_greedy   = NA_real_,
      cpu_exact    = NA_real_,
      error_msg    = e$message,
      stringsAsFactors = FALSE
    )
  })
}

# ------------------------------------------------------------------------------
# 5. Define the Scenario Grid (iter is the INNER parallel dimension)
# ------------------------------------------------------------------------------
scenario_grid <- expand.grid(
  N            = c(500, 1000, 2000, 5000),
  k            = c(1, 3, 5, 10, 15, 20),
  B_type       = c("20", "50", "100", "sqrt"),
  architecture = c("simple", "complex", "interaction",
                   "triple_interaction", "nonlinear_nuisance"),
  stringsAsFactors = FALSE
)

n_scenarios <- nrow(scenario_grid)
n_iters     <- sim_params$n_iters
total_rows  <- n_scenarios * n_iters

dir.create("../output/temp_03_03", recursive = TRUE, showWarnings = FALSE)

cat(sprintf(paste0(
  "Starting 03 Scaling Suite.\n",
  "  N    = {500, 1000, 2000, 5000}\n",
  "  k    = {1, 3, 5, 10, 20}\n",
  "  B    = {20, 50, 100, sqrt(N)}\n",
  "  Arch = {simple, complex, interaction, triple_interaction, nonlinear_nuisance}\n",
  "  Scenarios: %d | Iterations per scenario: %d | Total rows: %d\n\n"),
  n_scenarios, n_iters, total_rows))

# ------------------------------------------------------------------------------
# 6. Orchestrator Loop (Sequential scenarios, parallel iterations)
# ------------------------------------------------------------------------------
for (i in seq_len(n_scenarios)) {
  
  sc <- scenario_grid[i, ]
  chunk_file <- sprintf("../output/temp_03_03/03_chunk_%04d.rds", i)
  
  # Checkpoint: skip if already computed
  if (is_computed(chunk_file)) {
    if (i %% 20 == 0) cat(sprintf("[%04d/%04d] Skipping (cached): N=%d k=%d B=%s arch=%s\n",
                                  i, n_scenarios, sc$N, sc$k, sc$B_type, sc$architecture))
    next
  }
  
  cat(sprintf("[%04d/%04d] Computing: N=%d  k=%d  B=%s  arch=%s ... ",
              i, n_scenarios, sc$N, sc$k, sc$B_type, sc$architecture))
  
  # Resolve block count
  B <- if (sc$B_type == "sqrt") floor(sqrt(sc$N)) else as.numeric(sc$B_type)
  
  # Parallel inner loop over iterations
  scenario_results <- furrr::future_map_dfr(
    seq_len(n_iters), function(iter_id) {
      
      iter_seed <- sim_params$seed + (i - 1) * n_iters + iter_id
      
      run_scaling_iteration(
        iter_id      = iter_id,
        N            = sc$N,
        k            = sc$k,
        B            = B,
        architecture = sc$architecture,
        magnitude    = sim_params$magnitude,
        iter_seed    = iter_seed
      )
      
    }, .options = furrr_options(seed = TRUE)
  )
  
  # Attach scenario-level columns
  scenario_results$N            <- sc$N
  scenario_results$k            <- sc$k
  scenario_results$B_type       <- sc$B_type
  scenario_results$B_actual     <- B
  scenario_results$architecture <- sc$architecture
  
  safe_save_rds(scenario_results, chunk_file)
  cat("Done.\n")
}

# ------------------------------------------------------------------------------
# 7. Compile Final Artifact
# ------------------------------------------------------------------------------
cat("\nAll scenarios completed. Assembling final dataset...\n")
final_data <- compile_checkpoints(
  temp_dir          = "../output/temp_03_03",
  pattern           = "^03_chunk_.*\\.rds$",
  final_output_path = "../output/03_scaling_results_master.rds",
  clear_temp        = FALSE
)
cat("Script 03 execution finished successfully.\n")

# ------------------------------------------------------------------------------
# 8. Sanity Check & Reporting Validation
# ------------------------------------------------------------------------------
res <- final_data

cat("\n============================\n")
cat("  03 SANITY CHECK\n")
cat("============================\n")

cat("\n=== 1. BASIC INTEGRITY ===\n")
cat(sprintf("Rows: %d / %d (Missing: %d)\n",
            nrow(res), total_rows, total_rows - nrow(res)))
cat("Missing CPU (greedy):", sum(is.na(res$cpu_greedy)), "\n")
cat("Missing CPU (exact):",  sum(is.na(res$cpu_exact)), "\n")
cat("Missing detection rates:", sum(is.na(res$detection_rate)), "\n")

# Grid coverage: every (N, k, B_type, architecture) should have n_iters rows
coverage <- res %>%
  group_by(N, k, B_type, architecture) %>%
  summarise(n = n(), .groups = "drop")
incomplete <- coverage %>% filter(n < n_iters)
if (nrow(incomplete) > 0) {
  cat(sprintf("WARNING: %d scenarios have fewer than %d iterations:\n",
              nrow(incomplete), n_iters))
  print(incomplete)
} else {
  cat(sprintf("All %d scenarios have full %d iterations.\n", n_scenarios, n_iters))
}

cat("\n=== 2. STABILITY & ERROR RATES ===\n")
total_errors <- sum(!is.na(res$error_msg))
cat("Total algorithm crashes caught:", total_errors, "\n")

if (total_errors > 0) {
  cat("\nErrors by architecture:\n")
  print(table(res$architecture, Has_Error = !is.na(res$error_msg)))
  
  cat("\nErrors by (architecture, k):\n")
  print(table(res$architecture, res$k, !is.na(res$error_msg)))
  
  cat("\nSample error messages:\n")
  print(head(unique(na.omit(res$error_msg)), 3))
} else {
  cat("Perfect stability. No crashes.\n")
}

cat("\n=== 3. CONVERGENCE ===\n")
cat(sprintf("Greedy EVT converged: %.1f%%\n",
            mean(res$converged_greedy, na.rm = TRUE) * 100))
cat(sprintf("Exact  EVT converged: %.1f%%\n",
            mean(res$converged_exact, na.rm = TRUE) * 100))

conv_by_arch <- res %>%
  group_by(architecture) %>%
  summarise(
    conv_greedy = round(mean(converged_greedy, na.rm = TRUE) * 100, 1),
    conv_exact  = round(mean(converged_exact, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )
print(conv_by_arch)

conv_by_k <- res %>%
  group_by(k) %>%
  summarise(
    conv_greedy = round(mean(converged_greedy, na.rm = TRUE) * 100, 1),
    conv_exact  = round(mean(converged_exact, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )
cat("\nConvergence by k:\n")
print(conv_by_k)

conv_by_N <- res %>%
  group_by(N) %>%
  summarise(
    conv_greedy = round(mean(converged_greedy, na.rm = TRUE) * 100, 1),
    conv_exact  = round(mean(converged_exact, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )
cat("\nConvergence by N:\n")
print(conv_by_N)

cat("\n=== 4. DETECTION PREVIEW ===\n")
det_table <- res %>%
  filter(!is.na(detection_rate)) %>%
  group_by(architecture, k) %>%
  summarise(
    det_greedy = round(mean(detection_rate), 3),
    det_exact  = round(mean(detection_rate_exact, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(architecture, k)
print(det_table, n = 30)

# Warn on trivial results
all_zero <- all(det_table$det_greedy == 0) && all(det_table$det_exact == 0)
all_one  <- all(det_table$det_greedy == 1) && all(det_table$det_exact == 1)
if (all_zero) cat("WARNING: All detection rates are 0 — algorithm may be blind.\n")
if (all_one)  cat("WARNING: All detection rates are 1 — problem may be too easy.\n")

# Detection by N to see if larger samples help
det_by_N <- res %>%
  filter(!is.na(detection_rate)) %>%
  group_by(N, k) %>%
  summarise(
    det_greedy = round(mean(detection_rate), 3),
    det_exact  = round(mean(detection_rate_exact, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(N, k)
cat("\nDetection by (N, k):\n")
print(det_by_N, n = 25)

cat("\n=== 5. REPORTING SANITY ===\n")
cat("Greedy detection rates outside [0,1]:",
    sum(res$detection_rate < 0 | res$detection_rate > 1, na.rm = TRUE), "\n")
cat("Exact  detection rates outside [0,1]:",
    sum(res$detection_rate_exact < 0 | res$detection_rate_exact > 1, na.rm = TRUE), "\n")
cat("Negative CPU (greedy):", sum(res$cpu_greedy < 0, na.rm = TRUE), "\n")
cat("Negative CPU (exact):",  sum(res$cpu_exact < 0, na.rm = TRUE), "\n")
cat("Implausibly large CPU (>600s, greedy):",
    sum(res$cpu_greedy > 600, na.rm = TRUE), "\n")
cat("Implausibly large CPU (>600s, exact):",
    sum(res$cpu_exact > 600, na.rm = TRUE), "\n")

# B_actual consistency
res <- res %>%
  mutate(B_expected = ifelse(B_type == "sqrt", floor(sqrt(N)),
                             as.numeric(B_type)))
cat("B_actual mismatches:", sum(res$B_actual != res$B_expected, na.rm = TRUE), "\n")
res$B_expected <- NULL

# P-value sanity (converged rows should have p in [0,1])
p_bad_greedy <- res %>% filter(converged_greedy == TRUE & (p_greedy < 0 | p_greedy > 1))
p_bad_exact  <- res %>% filter(converged_exact  == TRUE & (p_exact  < 0 | p_exact  > 1))
cat("Greedy p-values outside [0,1] (converged):", nrow(p_bad_greedy), "\n")
cat("Exact  p-values outside [0,1] (converged):", nrow(p_bad_exact), "\n")

# NA in converged rows
na_greedy <- res %>% filter(converged_greedy == TRUE & is.na(p_greedy))
na_exact  <- res %>% filter(converged_exact  == TRUE & is.na(p_exact))
if (nrow(na_greedy) > 0) cat("WARNING:", nrow(na_greedy),
                             "converged greedy rows have NA p-values!\n")
if (nrow(na_exact) > 0)  cat("WARNING:", nrow(na_exact),
                             "converged exact rows have NA p-values!\n")

cat("\n=== 6. TIMING PROFILE ===\n")
timing_by_N <- res %>%
  filter(!is.na(cpu_greedy) & !is.na(cpu_exact)) %>%
  group_by(N) %>%
  summarise(
    med_greedy = round(median(cpu_greedy), 3),
    med_exact  = round(median(cpu_exact), 3),
    ratio      = round(median(cpu_exact) / max(median(cpu_greedy), 1e-6), 1),
    .groups = "drop"
  )
cat("By sample size (N):\n")
print(timing_by_N)

timing_by_k <- res %>%
  filter(!is.na(cpu_greedy) & !is.na(cpu_exact)) %>%
  group_by(k) %>%
  summarise(
    med_greedy = round(median(cpu_greedy), 3),
    med_exact  = round(median(cpu_exact), 3),
    ratio      = round(median(cpu_exact) / max(median(cpu_greedy), 1e-6), 1),
    .groups = "drop"
  )
cat("\nBy set size (k):\n")
print(timing_by_k)

timing_by_arch <- res %>%
  filter(!is.na(cpu_greedy) & !is.na(cpu_exact)) %>%
  group_by(architecture) %>%
  summarise(
    med_greedy = round(median(cpu_greedy), 3),
    med_exact  = round(median(cpu_exact), 3),
    ratio      = round(median(cpu_exact) / max(median(cpu_greedy), 1e-6), 1),
    .groups = "drop"
  )
cat("\nBy architecture:\n")
print(timing_by_arch)

# Cross-tabulation: N x k (the most informative scaling view)
timing_Nk <- res %>%
  filter(!is.na(cpu_exact)) %>%
  group_by(N, k) %>%
  summarise(
    med_exact = round(median(cpu_exact), 3),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = k, values_from = med_exact,
                     names_prefix = "k=")
cat("\nExact median CPU (sec) — N x k:\n")
print(timing_Nk)

cat("\nSanity check complete.\n")