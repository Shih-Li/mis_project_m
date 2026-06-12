# ==============================================================================
# File: scripts/03_alg_comp_v2.R
# Purpose: Revised algorithmic comparison with 2x2 detection factorial
#          ({greedy, exact} x {no refinement, with refinement}) and 3 EVT
#          configurations ({greedy+constrained, exact+constrained, exact+robust}).
#          Adds plm_confounded and plm_nonlinear architectures.
#          Uses new Phase 0 functions and metrics_influence.R.
#
# Outputs:
#   ../output/temp_03v2/03v2_chunk_*.rds  (per-scenario checkpoints)
#   ../output/03v2_scaling_results_master.rds  (compiled main results)
#   ../output/03v2_nestedness_traces.rds  (nestedness sub-analysis)
# ==============================================================================

# 1. Load Dependencies
library(dplyr)
library(purrr)
library(evd)
library(future)
library(furrr)

source("../R/helpers_local.R")
source("../R/03_scaling_dgp.R")       # includes plm_confounded, plm_nonlinear
source("../R/utils_checkpoint.R")
source("../R/fast_sens_topk.R")        # includes fast_sens_topk_diag
source("../R/dinkelbach_topk.R")       # includes dinkelbach_topk_refined
source("../R/exact_dfb_bmx.R")         # includes exact_dfb_bmx_diag
source("../R/evt_iter.R")
source("../R/evt_iter_dm.R")           # includes evt_iter_dm_v2
source("../R/metrics_influence.R")
source("../R/sim_scaling_engine_v2.R")

# Parallel setup
slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK")
if (slurm_cpus != "") {
  num_workers <- as.numeric(slurm_cpus)
  cat(sprintf("SLURM environment, using %d workers.\n", num_workers))
} else {
  num_workers <- max(1, parallel::detectCores() - 2)
  cat(sprintf("Local environment, using %d workers.\n", num_workers))
}
plan(multisession, workers = num_workers)

# ------------------------------------------------------------------------------
# 2. Global Configuration
# ------------------------------------------------------------------------------
sim_params <- list(
  n_iters   = 100,
  magnitude = 10,
  seed      = 20260610
)

set.seed(sim_params$seed)

MIN_OBS_PER_BLOCK <- 30

resolve_block_count <- function(B_type, N, k) {
  B_raw <- if (B_type == "sqrt") floor(sqrt(N)) else as.numeric(B_type)
  B_max <- floor((N - k) / MIN_OBS_PER_BLOCK)
  B_max <- max(B_max, 3L)
  min(B_raw, B_max)
}

# ------------------------------------------------------------------------------
# 3. Scenario Grid
# ------------------------------------------------------------------------------

# --- Base grid (existing architectures) ---
grid_base <- expand.grid(
  N            = c(500, 1000, 2000, 5000),
  k            = c(1, 3, 5, 10, 15, 20),
  B_type       = c("20", "50", "100", "sqrt"),
  architecture = c("simple", "complex", "interaction",
                   "triple_interaction", "nonlinear_nuisance",
                   "sparse_binary_interaction",
                   "polynomial_interaction",
                   "plm_confounded", "plm_nonlinear"),
  stringsAsFactors = FALSE
)

# --- high_k_interaction ---
grid_high_k <- expand.grid(
  N            = c(500, 1000, 2000, 5000),
  k            = c(1, 3, 5, 10, 15, 20, 30, 50),
  B_type       = c("20", "50", "100", "sqrt"),
  architecture = "high_k_interaction",
  stringsAsFactors = FALSE
)

# --- Rho sweep for collinear_interaction ---
grid_collinear <- expand.grid(
  N            = c(500, 1000, 2000, 5000),
  k            = c(1, 3, 5, 10, 15, 20),
  B_type       = c("20", "50", "100", "sqrt"),
  architecture = "collinear_interaction",
  rho          = c(0.5, 0.7, 0.85, 0.95),
  stringsAsFactors = FALSE
)

# Combine
grid_base$rho    <- NA_real_
grid_high_k$rho  <- NA_real_
scenario_grid <- rbind(grid_base, grid_high_k, grid_collinear)
scenario_grid <- scenario_grid[scenario_grid$k / scenario_grid$N <= 0.05, ]
scenario_grid <- scenario_grid[order(scenario_grid$architecture,
                                     scenario_grid$N,
                                     scenario_grid$k), ]
rownames(scenario_grid) <- NULL

n_scenarios <- nrow(scenario_grid)
n_iters     <- sim_params$n_iters
total_rows  <- n_scenarios * n_iters

dir.create("../output/temp_03v2", recursive = TRUE, showWarnings = FALSE)

cat(sprintf(paste0(
  "Starting 03v2 Scaling Suite (2x2 Detection Factorial).\n",
  "  Architectures: %d (incl. plm_confounded, plm_nonlinear)\n",
  "  Scenarios: %d | Iterations: %d | Total rows: %d\n\n"),
  length(unique(scenario_grid$architecture)),
  n_scenarios, n_iters, total_rows))

# ------------------------------------------------------------------------------
# 4. Orchestrator Loop
# ------------------------------------------------------------------------------
for (i in seq_len(n_scenarios)) {
  
  sc <- scenario_grid[i, ]
  chunk_file <- sprintf("../output/temp_03v2/03v2_chunk_%04d.rds", i)
  
  if (is_computed(chunk_file)) {
    if (i %% 50 == 0) cat(sprintf("[%04d/%04d] Skipping (cached)\n",
                                  i, n_scenarios))
    next
  }
  
  cat(sprintf("[%04d/%04d] N=%d k=%d B=%s arch=%s ... ",
              i, n_scenarios, sc$N, sc$k, sc$B_type, sc$architecture))
  
  B <- resolve_block_count(sc$B_type, sc$N, sc$k)
  B_raw <- if (sc$B_type == "sqrt") floor(sqrt(sc$N)) else as.numeric(sc$B_type)
  B_was_capped <- (B < B_raw)
  if (B_was_capped) cat(sprintf("[B: %d->%d] ", B_raw, B))
  
  scenario_results <- furrr::future_map_dfr(
    seq_len(n_iters), function(iter_id) {
      
      iter_seed <- sim_params$seed + (i - 1) * n_iters + iter_id
      
      run_scaling_iteration_v2(
        iter_id      = iter_id,
        N            = sc$N,
        k            = sc$k,
        B            = B,
        architecture = sc$architecture,
        magnitude    = sim_params$magnitude,
        iter_seed    = iter_seed,
        rho          = sc$rho
      )
      
    }, .options = furrr_options(seed = TRUE)
  )
  
  # Attach scenario-level metadata
  scenario_results$N            <- sc$N
  scenario_results$k            <- sc$k
  scenario_results$B_type       <- sc$B_type
  scenario_results$B_actual     <- B
  scenario_results$B_raw        <- B_raw
  scenario_results$B_capped     <- B_was_capped
  scenario_results$architecture <- sc$architecture
  scenario_results$rho          <- sc$rho
  
  safe_save_rds(scenario_results, chunk_file)
  cat("Done.\n")
}

# ------------------------------------------------------------------------------
# 5. Compile Final Artifact
# ------------------------------------------------------------------------------
cat("\nAssembling final dataset...\n")
final_data <- compile_checkpoints(
  temp_dir          = "../output/temp_03v2",
  pattern           = "^03v2_chunk_.*\\.rds$",
  final_output_path = "../output/03v2_scaling_results_master.rds",
  clear_temp        = FALSE
)
cat("Compiled successfully.\n")

# ------------------------------------------------------------------------------
# 6. Nestedness Sub-Analysis (representative subset)
# ------------------------------------------------------------------------------
cat("\nRunning nestedness traces...\n")

nest_archs <- c("simple", "interaction", "collinear_interaction",
                "plm_nonlinear")
nest_N     <- 2000
nest_k_max <- 20

set.seed(sim_params$seed + 999)

nest_results <- list()
counter <- 0L

for (arch in nest_archs) {
  dgp <- generate_scaling_dgp(N = nest_N, architecture = arch)
  dgp_inj <- inject_safe_outliers_v2(dgp, k = nest_k_max,
                                     magnitude = sim_params$magnitude)
  mod <- lm(dgp_inj$formula, data = dgp_inj$data)
  tpos <- dgp_inj$target_pos
  
  for (method in c("greedy", "dinkelbach", "dinkelbach_refined")) {
    counter <- counter + 1L
    cat(sprintf("  Nestedness: arch=%s method=%s ... ", arch, method))
    
    trace <- tryCatch(
      nestedness_trace(mod, pos = tpos, k_max = nest_k_max,
                       method = method, sign = 1L),
      error = function(e) {
        warning(sprintf("Nestedness failed: %s / %s: %s",
                        arch, method, e$message))
        NULL
      }
    )
    
    if (!is.null(trace)) {
      trace$architecture <- arch
      trace$method       <- method
      nest_results[[counter]] <- trace
    }
    cat("Done.\n")
  }
}

nest_df <- do.call(rbind, nest_results)
safe_save_rds(nest_df, "../output/03v2_nestedness_traces.rds")
cat(sprintf("Nestedness traces saved (%d rows).\n", nrow(nest_df)))

# ------------------------------------------------------------------------------
# 7. Sanity Check
# ------------------------------------------------------------------------------
res <- final_data

cat("\n============================\n")
cat("  03v2 SANITY CHECK\n")
cat("============================\n")

cat("\n=== 1. BASIC INTEGRITY ===\n")
cat(sprintf("Rows: %d / %d expected\n", nrow(res), total_rows))
cat("Errors caught:", sum(!is.na(res$error_msg)), "\n")

cat("\n=== 2. DETECTION FACTORIAL ===\n")
det_table <- res %>%
  filter(!is.na(det_greedy_noref)) %>%
  group_by(architecture, k) %>%
  summarise(
    greedy_noref = round(mean(det_greedy_noref), 3),
    greedy_ref   = round(mean(det_greedy_ref), 3),
    exact_noref  = round(mean(det_exact_noref), 3),
    exact_ref    = round(mean(det_exact_ref), 3),
    .groups = "drop"
  ) %>%
  arrange(architecture, k)
cat("Mean detection overlap by (architecture, k):\n")
print(det_table, n = 40)

cat("\n=== 3. INFLUENCE RATIO ===\n")
ir_table <- res %>%
  filter(!is.na(ir_greedy_noref)) %>%
  group_by(architecture) %>%
  summarise(
    greedy_noref = round(mean(ir_greedy_noref, na.rm = TRUE), 3),
    greedy_ref   = round(mean(ir_greedy_ref, na.rm = TRUE), 3),
    exact_noref  = round(mean(ir_exact_noref, na.rm = TRUE), 3),
    exact_ref    = round(mean(ir_exact_ref, na.rm = TRUE), 3),
    .groups = "drop"
  )
cat("Mean influence ratio by architecture:\n")
print(ir_table)

cat("\n=== 4. EVT CONVERGENCE (3 configurations) ===\n")
evt_table <- res %>%
  group_by(architecture) %>%
  summarise(
    conv_greedy_con = round(mean(conv_greedy_con, na.rm = TRUE) * 100, 1),
    conv_exact_con  = round(mean(conv_exact_con, na.rm = TRUE) * 100, 1),
    conv_exact_rob  = round(mean(conv_exact_rob, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )
cat("Convergence rates (%):\n")
print(evt_table)

cat("\n=== 5. EVT SHAPE PARAMETERS ===\n")
shape_table <- res %>%
  filter(conv_exact_con == TRUE) %>%
  group_by(architecture) %>%
  summarise(
    shape_greedy_con = round(mean(shape_greedy_con, na.rm = TRUE), 3),
    shape_exact_con  = round(mean(shape_exact_con, na.rm = TRUE), 3),
    shape_exact_rob  = round(mean(shape_exact_rob, na.rm = TRUE), 3),
    .groups = "drop"
  )
cat("Mean shape (xi) by architecture:\n")
print(shape_table)

cat("\n=== 6. TIMING ===\n")
timing <- res %>%
  filter(!is.na(cpu_det_greedy_noref)) %>%
  group_by(N) %>%
  summarise(
    det_gn = round(median(cpu_det_greedy_noref), 4),
    det_gr = round(median(cpu_det_greedy_ref), 4),
    det_en = round(median(cpu_det_exact_noref), 4),
    det_er = round(median(cpu_det_exact_ref), 4),
    evt_gc = round(median(cpu_evt_greedy_con, na.rm = TRUE), 4),
    evt_ec = round(median(cpu_evt_exact_con, na.rm = TRUE), 4),
    evt_er = round(median(cpu_evt_exact_rob, na.rm = TRUE), 4),
    .groups = "drop"
  )
cat("Median CPU (sec) by N:\n")
print(timing)

cat("\n=== 7. PLM ARCHITECTURES ===\n")
plm_table <- res %>%
  filter(architecture %in% c("plm_confounded", "plm_nonlinear")) %>%
  group_by(architecture, k) %>%
  summarise(
    det_exact_ref = round(mean(det_exact_ref), 3),
    ir_exact_ref  = round(mean(ir_exact_ref, na.rm = TRUE), 3),
    conv_exact_con = round(mean(conv_exact_con, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )
if (nrow(plm_table) > 0) {
  cat("PLM detection + EVT:\n")
  print(plm_table, n = 30)
}

cat("\n=== 8. KEY HYPOTHESIS CHECK ===\n")
h1 <- res %>%
  filter(!is.na(det_greedy_ref)) %>%
  summarise(
    greedy_ref_mean = round(mean(det_greedy_ref), 3),
    exact_ref_mean  = round(mean(det_exact_ref), 3),
    exact_noref_mean = round(mean(det_exact_noref), 3),
    greedy_noref_mean = round(mean(det_greedy_noref), 3)
  )
cat("H1 (refinement is the key driver):\n")
cat(sprintf("  With refinement:    greedy=%.3f  exact=%.3f  (should be similar)\n",
            h1$greedy_ref_mean, h1$exact_ref_mean))
cat(sprintf("  Without refinement: greedy=%.3f  exact=%.3f  (both worse)\n",
            h1$greedy_noref_mean, h1$exact_noref_mean))

cat("\nScript 03v2 complete.\n")