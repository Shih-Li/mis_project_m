# ==============================================================================
# File: scripts/05_misspecification.R
# Purpose:
#   Controlled model-misspecification simulation.
#
# Scientific design:
#
#   N
#   x distribution
#   error distribution
#   misspecification type
#   misspecification severity
#
# Each outer design cell is:
#
#   N × X distribution × error distribution × scenario × severity
#
# All fixed deletion fractions k/n are evaluated INSIDE the same
# Monte Carlo draw:
#
#   0.01
#   0.025
#   0.05
#   0.10
#
# Diagnostics:
#
#   MIS
#   DFBETAS
#   Cook's distance
#   Leverage
#
# Estimators:
#
#   Full OLS
#   Full MM
#   MIS -> OLS
#   MIS -> MM
#   Classical diagnostic -> OLS
#
# Execution:
#
#   - Monte Carlo iterations parallelised across CPU workers
#   - optional Rcpp Dinkelbach kernel
#   - checkpoint saved once per completed outer design cell
#   - common random numbers across severity levels
#
# IMPORTANT:
#
#   Severity levels for the same
#
#       N × X × error × scenario × iteration
#
#   use the same underlying RNG seed.
#
#   Therefore changes across severity are driven primarily by the
#   misspecification-strength parameter rather than by a completely
#   different random sample.
# ==============================================================================


# ==============================================================================
# 1. Packages
# ==============================================================================

library(dplyr)
library(future)
library(furrr)


# ==============================================================================
# 2. Project functions
# ==============================================================================

source("../R/helpers_local.R")

source("../R/diagnostics_classical.R")
source("../R/dinkelbach_topk.R")

source("../R/dgp_misspec_factory.R")
source("../R/mis_sensitivity.R")
source("../R/misspec_metrics.R")
source("../R/sim_misspec_engine.R")

# Existing checkpoint infrastructure from Scripts 02/04
source("../R/utils_checkpoint.R")


# ==============================================================================
# 3. Simulation configuration
# ==============================================================================

sim_params <- list(
  
  # --------------------------------------------------------------------------
  # Monte Carlo
  # --------------------------------------------------------------------------
  
  n_iters = 100L,
  
  seed = 20260828L,
  
  
  # --------------------------------------------------------------------------
  # C++ acceleration
  # --------------------------------------------------------------------------
  
  use_cpp = TRUE,
  
  
  # --------------------------------------------------------------------------
  # Parallel workers
  #
  # Leave one logical core free and cap at 12 workers initially.
  # --------------------------------------------------------------------------
  
  workers = max(
    1L,
    min(
      future::availableCores() - 1L,
      12L
    )
  ),
  
  
  # --------------------------------------------------------------------------
  # Distribution parameters
  # --------------------------------------------------------------------------
  
  # Main GPD experiment:
  #
  # shape < 0.5 gives finite variance
  gpd_shape = 0.25,
  
  pareto_shape = 3,
  
  skew_t_df = 5,
  
  
  # --------------------------------------------------------------------------
  # OVB
  # --------------------------------------------------------------------------
  
  rho_xz = 0.50,
  
  
  # --------------------------------------------------------------------------
  # Heterogeneous slope
  # --------------------------------------------------------------------------
  
  hetero_group_prop = 0.25,
  
  
  # --------------------------------------------------------------------------
  # Structural break
  # --------------------------------------------------------------------------
  
  break_fraction = 0.75,
  
  
  # --------------------------------------------------------------------------
  # Threshold model
  # --------------------------------------------------------------------------
  
  threshold_quantile = 0.75
)


# ==============================================================================
# 4. Experimental grid
# ==============================================================================

n_grid <- c(
  500L,
  1000L,
  2500L,
  5000L
)


# ------------------------------------------------------------------------------
# Keep same X environments as previous robustness simulation
# ------------------------------------------------------------------------------

x_grid <- c(
  "normal",
  "mixed_normal",
  "contaminated"
)


# ------------------------------------------------------------------------------
# Keep same eight error environments
# ------------------------------------------------------------------------------

error_grid <- c(
  "normal",
  "mixed_normal",
  "skewed_t",
  "golm",
  "beta_logistic",
  "gpd",
  "contaminated",
  "pareto"
)


# ------------------------------------------------------------------------------
# Structural/model problems investigated
# ------------------------------------------------------------------------------

scenario_grid <- c(
  
  "correct",
  
  "ovb",
  
  "nonlinear",
  
  "heterogeneous",
  
  "structural_break",
  
  "threshold",
  
  "missing_interaction",
  
  # Different-mechanism / negative controls
  "endogeneity",
  
  "heteroskedastic"
)


# ------------------------------------------------------------------------------
# Severity level
#
# Additive structural misspecification:
#
#   level 0 -> lambda = 0
#   level 1 -> lambda = 0.25
#   level 2 -> lambda = 0.50
#   level 3 -> lambda = 1.00
#   level 4 -> lambda = 2.00
#
# Endogeneity and heteroskedasticity use scenario-specific mappings
# implemented in dgp_misspec_factory.R.
# ------------------------------------------------------------------------------

severity_level_grid <- 0:4


# ------------------------------------------------------------------------------
# Fixed deletion budget -- NOT oracle k
#
# IMPORTANT:
#
# These are evaluated inside each Monte Carlo dataset.
# They are NOT separate outer design cells.
# ------------------------------------------------------------------------------

k_fraction_grid <- c(
  0.01,
  0.025,
  0.05,
  0.10
)


# ==============================================================================
# 5. Construct design cells
# ==============================================================================

design_grid <- expand.grid(
  
  n = n_grid,
  
  x_type = x_grid,
  
  error_type = error_grid,
  
  scenario = scenario_grid,
  
  severity_level = severity_level_grid,
  
  stringsAsFactors = FALSE
  
) |>
  
  # --------------------------------------------------------------------------
# Correct model has only severity = 0
# --------------------------------------------------------------------------

filter(
  !(
    scenario == "correct" &
      severity_level != 0L
  )
) |>
  
  # --------------------------------------------------------------------------
# Stable ordering
# --------------------------------------------------------------------------

arrange(
  n,
  x_type,
  error_type,
  scenario,
  severity_level
) |>
  
  # --------------------------------------------------------------------------
# COMMON RANDOM NUMBER GROUP
#
# Severity is intentionally NOT included in this grouping.
#
# Hence:
#
# N = 500
# X = normal
# error = normal
# scenario = nonlinear
#
# gets ONE seed_group_id shared by:
#
# severity 0
# severity 1
# severity 2
# severity 3
# severity 4
# --------------------------------------------------------------------------

group_by(
  n,
  x_type,
  error_type,
  scenario
) |>
  
  mutate(
    seed_group_id =
      dplyr::cur_group_id()
  ) |>
  
  ungroup() |>
  
  # --------------------------------------------------------------------------
# Unique checkpoint cell ID
#
# Unlike seed_group_id, cell_id DOES distinguish severity.
# --------------------------------------------------------------------------

mutate(
  cell_id =
    row_number()
)


message(
  "Script 05 design cells: ",
  nrow(design_grid)
)


message(
  "Common-RNG seed groups: ",
  dplyr::n_distinct(
    design_grid$seed_group_id
  )
)


# ==============================================================================
# 6. Output / checkpoint directories
# ==============================================================================

output_dir <- "../output"


# IMPORTANT:
# Keep existing checkpoint directory name unchanged.

checkpoint_dir <- file.path(
  output_dir,
  "temp",
  "05_misspecification"
)


dir.create(
  checkpoint_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Legacy single-file raw-results path.
# Intentionally not written because the complete formal result is too large
# to compile safely in memory. Raw results remain in checkpoint RDS files.

final_results_path <- file.path(
  output_dir,
  "05_misspecification_results.rds"
)


summary_path <- file.path(
  output_dir,
  "05_misspecification_summary.rds"
)


power_path <- file.path(
  output_dir,
  "05_detection_power.rds"
)


boundary_path <- file.path(
  output_dir,
  "05_detection_boundary.rds"
)


# ==============================================================================
# 7. Parallel workers
# ==============================================================================
#
# Use an explicit PSOCK cluster.
#
# Why?
#
# multisession/cluster workers are separate R processes.
# A sourceCpp() call in the main R session does NOT guarantee that the
# compiled function exists inside every worker.
#
# Therefore:
#
# create workers
#      ->
# compile C++ once inside every worker
#      ->
# run simulation
# ==============================================================================

n_workers <- sim_params$workers


message(
  "Starting ",
  n_workers,
  " parallel workers."
)


cl <- parallel::makePSOCKcluster(
  n_workers
)


# ==============================================================================
# 8. Optional C++ MIS kernel
# ==============================================================================

cpp_ok <- FALSE


if (
  isTRUE(sim_params$use_cpp) &&
  requireNamespace(
    "Rcpp",
    quietly = TRUE
  ) &&
  file.exists(
    "../src/dinkelbach_topk_cpp.cpp"
  )
) {
  
  cpp_path <- normalizePath(
    "../src/dinkelbach_topk_cpp.cpp",
    mustWork = TRUE
  )
  
  
  # --------------------------------------------------------------------------
  # Make C++ source path available to workers
  # --------------------------------------------------------------------------
  
  parallel::clusterExport(
    cl,
    "cpp_path",
    envir = environment()
  )
  
  
  # --------------------------------------------------------------------------
  # Compile/load once per worker
  # --------------------------------------------------------------------------
  
  cpp_status <- parallel::clusterCall(
    cl,
    function() {
      
      if (
        !requireNamespace(
          "Rcpp",
          quietly = TRUE
        )
      ) {
        
        return(
          FALSE
        )
      }
      
      
      tryCatch(
        {
          
          Rcpp::sourceCpp(
            cpp_path,
            rebuild = FALSE,
            showOutput = FALSE,
            verbose = FALSE
          )
          
          
          exists(
            "dinkelbach_topk_cpp",
            mode = "function"
          )
          
        },
        error = function(e) {
          
          FALSE
        }
      )
    }
  )
  
  
  cpp_ok <- all(
    unlist(
      cpp_status
    )
  )
}


message(
  "C++ MIS kernel: ",
  if (cpp_ok) {
    "enabled on all workers"
  } else {
    "disabled; using pure-R Dinkelbach"
  }
)


# ------------------------------------------------------------------------------
# Use explicit worker cluster
# ------------------------------------------------------------------------------

future::plan(
  future::cluster,
  workers = cl
)


message(
  "Active workers: ",
  future::nbrOfWorkers()
)


# ==============================================================================
# 9. Run simulation cell by cell
# ==============================================================================
#
# CHECKPOINT LOGIC
#
# Each OUTER design cell is:
#
#   N
#   × X distribution
#   × error distribution
#   × scenario
#   × severity
#
# It does NOT include k/n.
#
# Inside one Monte Carlo draw, all:
#
#   1%
#   2.5%
#   5%
#   10%
#
# deletion fractions are evaluated using the same generated dataset.
#
#
# COMMON RANDOM NUMBERS
#
# For a given:
#
#   N
#   X distribution
#   error distribution
#   scenario
#   iteration
#
# all severity levels receive the same seed_base.
#
#
# CHECKPOINT:
#
# Once all repetitions for one outer design cell finish:
#
#   safe_save_rds()
#
# creates:
#
#   ../output/temp/05_misspecification/
#
#       05_cell_000001.rds
#       05_cell_000002.rds
#       ...
#
# If Script 05 is restarted, is_computed() skips completed cells.
# ==============================================================================


tryCatch({
  
  for (
    i in seq_len(
      nrow(design_grid)
    )
  ) {
    
    
    cell <-
      design_grid[
        i,
        ,
        drop = FALSE
      ]
    
    
    checkpoint_file <- file.path(
      checkpoint_dir,
      sprintf(
        "05_cell_%06d.rds",
        cell$cell_id
      )
    )
    
    
    # =========================================================================
    # Resume support
    # =========================================================================
    
    if (
      is_computed(
        checkpoint_file
      )
    ) {
      
      message(
        sprintf(
          "[%d/%d] Cell %d already completed -- skipping.",
          i,
          nrow(design_grid),
          cell$cell_id
        )
      )
      
      
      next
    }
    
    
    # =========================================================================
    # Progress message
    # =========================================================================
    
    message(
      sprintf(
        paste0(
          "[%d/%d] ",
          "N=%d | X=%s | error=%s | ",
          "scenario=%s | severity=%d | ",
          "seed-group=%d"
        ),
        
        i,
        
        nrow(design_grid),
        
        cell$n,
        
        cell$x_type,
        
        cell$error_type,
        
        cell$scenario,
        
        cell$severity_level,
        
        cell$seed_group_id
      )
    )
    
    
    # =========================================================================
    # Common random-number seed
    #
    # CRITICAL:
    #
    # seed_group_id excludes severity.
    #
    # Hence all severity levels for the same environment receive exactly
    # the same seed_base.
    #
    # run_misspec_iter() then uses:
    #
    #   set.seed(seed_base + iter)
    #
    # so the same iteration is matched across severity.
    # =========================================================================
    
    seed_base <- as.integer(
      sim_params$seed +
        100000L *
        as.integer(
          cell$seed_group_id
        )
    )
    
    
    # =========================================================================
    # Monte Carlo iterations in parallel
    # =========================================================================
    
    cell_results <- furrr::future_map_dfr(
      
      seq_len(
        sim_params$n_iters
      ),
      
      function(iter) {
        
        run_misspec_iter(
          
          # -------------------------------------------------------------------
          # Monte Carlo
          # -------------------------------------------------------------------
          
          iter =
            iter,
          
          
          # -------------------------------------------------------------------
          # Design
          # -------------------------------------------------------------------
          
          n =
            cell$n,
          
          x_type =
            cell$x_type,
          
          error_type =
            cell$error_type,
          
          scenario =
            cell$scenario,
          
          severity_level =
            cell$severity_level,
          
          
          # -------------------------------------------------------------------
          # All deletion fractions evaluated inside this same draw
          # -------------------------------------------------------------------
          
          k_fractions =
            k_fraction_grid,
          
          
          # -------------------------------------------------------------------
          # Common RNG seed
          # -------------------------------------------------------------------
          
          seed_base =
            seed_base,
          
          
          # -------------------------------------------------------------------
          # C++ acceleration
          # -------------------------------------------------------------------
          
          use_cpp =
            cpp_ok,
          
          
          # -------------------------------------------------------------------
          # DGP controls
          # -------------------------------------------------------------------
          
          rho_xz =
            sim_params$rho_xz,
          
          hetero_group_prop =
            sim_params$hetero_group_prop,
          
          break_fraction =
            sim_params$break_fraction,
          
          threshold_quantile =
            sim_params$threshold_quantile,
          
          
          # -------------------------------------------------------------------
          # Distribution controls
          # -------------------------------------------------------------------
          
          gpd_shape =
            sim_params$gpd_shape,
          
          pareto_shape =
            sim_params$pareto_shape,
          
          skew_t_df =
            sim_params$skew_t_df
        )
      },
      
      
      .options =
        furrr::furrr_options(
          
          # furrr itself remains reproducible.
          #
          # run_misspec_iter() additionally sets the explicit matched
          # seed_base + iter seed.
          seed = TRUE,
          
          scheduling = 2
        )
    )
    
    
    # =========================================================================
    # Explicit outer-design identifiers
    # =========================================================================
    
    cell_results$cell_id <-
      cell$cell_id
    
    
    cell_results$seed_group_id <-
      cell$seed_group_id
    
    
    # =========================================================================
    # Atomic checkpoint save
    # =========================================================================
    
    safe_save_rds(
      cell_results,
      checkpoint_file
    )
    
    
    message(
      sprintf(
        "Cell %d saved.",
        cell$cell_id
      )
    )
  }
  
  
  # ============================================================================
  # 10. Memory-safe checkpoint post-processing
  #
  # IMPORTANT:
  #
  # The formal experiment contains approximately 19 million raw result rows.
  # Do NOT compile all checkpoint rows into one in-memory data frame.
  #
  # The 3936 checkpoint RDS files are the raw formal output.
  #
  # This section reads one checkpoint at a time and constructs only the much
  # smaller:
  #
  #   05_misspecification_summary.rds
  #   05_detection_power.rds
  #   05_detection_boundary.rds
  #
  # The raw checkpoint files remain unchanged.
  # ============================================================================
  
  
  checkpoint_files <- sort(
    list.files(
      checkpoint_dir,
      pattern = "^05_cell_[0-9]+\\.rds$",
      full.names = TRUE
    )
  )
  
  
  if (
    length(checkpoint_files) !=
    nrow(design_grid)
  ) {
    
    stop(
      "Checkpoint count mismatch: expected ",
      nrow(design_grid),
      ", found ",
      length(checkpoint_files),
      "."
    )
  }
  
  
  message(
    "Post-processing ",
    length(checkpoint_files),
    " completed checkpoints."
  )
  
  
  # ============================================================================
  # 11. Safe summary helper functions
  # ============================================================================
  
  mean_or_na05 <- function(x) {
    
    x <-
      x[
        is.finite(x)
      ]
    
    
    if (
      length(x) == 0L
    ) {
      
      return(
        NA_real_
      )
    }
    
    
    mean(
      x
    )
  }
  
  
  rmse_or_na05 <- function(x) {
    
    x <-
      x[
        is.finite(x)
      ]
    
    
    if (
      length(x) == 0L
    ) {
      
      return(
        NA_real_
      )
    }
    
    
    sqrt(
      mean(
        x^2
      )
    )
  }
  
  
  # ============================================================================
  # 12. Detailed summary table — checkpoint by checkpoint
  #
  # Every grouping combination below lives entirely inside one checkpoint,
  # because a checkpoint already fixes:
  #
  #   n
  #   x_type
  #   error_type
  #   scenario
  #   severity_level
  #
  # Therefore each checkpoint can be summarised independently and the small
  # summaries can safely be combined afterward.
  # ============================================================================
  
  
  summary_parts <-
    vector(
      "list",
      length(checkpoint_files)
    )
  
  
  checkpoint_row_counts <-
    numeric(
      length(checkpoint_files)
    )
  
  
  for (
    j in seq_along(checkpoint_files)
  ) {
    
    x <-
      readRDS(
        checkpoint_files[j]
      )
    
    
    checkpoint_row_counts[j] <-
      nrow(x)
    
    
    summary_parts[[j]] <-
      x |>
      
      group_by(
        
        n,
        
        x_type,
        
        error_type,
        
        scenario,
        
        severity_level,
        
        severity_target,
        
        k_fraction,
        
        model_state,
        
        diagnostic,
        
        estimator
        
      ) |>
      
      summarise(
        
        n_records =
          dplyr::n(),
        
        
        true_beta =
          mean_or_na05(
            true_beta
          ),
        
        
        mean_coef =
          mean_or_na05(
            coef
          ),
        
        mean_se =
          mean_or_na05(
            se
          ),
        
        mean_bias =
          mean_or_na05(
            bias
          ),
        
        mean_abs_bias =
          mean_or_na05(
            abs_bias
          ),
        
        rmse =
          rmse_or_na05(
            bias
          ),
        
        coverage =
          mean_or_na05(
            coverage
          ),
        
        
        mean_delta_beta =
          mean_or_na05(
            delta_beta
          ),
        
        mean_standardized_delta =
          mean_or_na05(
            standardized_delta
          ),
        
        
        mean_signed_delta_beta =
          mean_or_na05(
            signed_delta_beta
          ),
        
        mean_signed_standardized_delta =
          mean_or_na05(
            signed_standardized_delta
          ),
        
        mean_relative_delta_beta =
          mean_or_na05(
            relative_delta_beta
          ),
        
        
        mean_outlier_rate =
          mean_or_na05(
            outlier_rate
          ),
        
        mean_n_abnormal =
          mean_or_na05(
            n_abnormal
          ),
        
        mean_residual_mad =
          mean_or_na05(
            residual_mad
          ),
        
        mean_q95_abs_std_resid =
          mean_or_na05(
            q95_abs_std_resid
          ),
        
        mean_max_abs_std_resid =
          mean_or_na05(
            max_abs_std_resid
          ),
        
        
        mean_max_leverage =
          mean_or_na05(
            max_leverage
          ),
        
        mean_condition_number =
          mean_or_na05(
            condition_number
          ),
        
        
        mean_affected_n =
          mean_or_na05(
            affected_n
          ),
        
        mean_affected_fraction =
          mean_or_na05(
            affected_fraction
          ),
        
        mean_overlap_n =
          mean_or_na05(
            overlap_n
          ),
        
        mean_overlap_recall =
          mean_or_na05(
            overlap_recall
          ),
        
        mean_overlap_precision =
          mean_or_na05(
            overlap_precision
          ),
        
        mean_selection_lift =
          mean_or_na05(
            selection_lift
          ),
        
        
        mean_selected_size =
          mean_or_na05(
            selected_size
          ),
        
        
        mean_runtime =
          mean_or_na05(
            runtime
          ),
        
        
        .groups =
          "drop"
      )
    
    
    rm(x)
    
    
    if (
      j %% 250L ==
      0L
    ) {
      
      gc(
        verbose = FALSE
      )
    }
  }
  
  
  summary_table <-
    dplyr::bind_rows(
      summary_parts
    )
  
  
  rm(
    summary_parts
  )
  
  
  gc(
    verbose = FALSE
  )
  
  
  safe_save_rds(
    summary_table,
    summary_path
  )
  
  
  message(
    "Summary table saved: ",
    summary_path
  )
  
  
  # ============================================================================
  # 13. Detection power — memory-safe two-pass calculation
  #
  # Existing definition:
  #
  #   1. severity_level = 0 determines the 95% null cutoff
  #   2. each severity is compared with its matching null cutoff
  #
  # Matching dimensions:
  #
  #   n
  #   x_type
  #   error_type
  #   scenario
  #   k_fraction
  #   diagnostic
  # ============================================================================
  
  
  detection_group_vars <- c(
    "n",
    "x_type",
    "error_type",
    "scenario",
    "k_fraction",
    "diagnostic"
  )
  
  
  # ----------------------------------------------------------------------------
  # PASS 1:
  # Calculate severity-0 null cutoffs
  # ----------------------------------------------------------------------------
  
  null_parts <-
    vector(
      "list",
      length(checkpoint_files)
    )
  
  
  for (
    j in seq_along(checkpoint_files)
  ) {
    
    x <-
      readRDS(
        checkpoint_files[j]
      )
    
    
    if (
      unique(
        x$severity_level
      )[1L] ==
      0L
    ) {
      
      z <-
        x |>
        
        filter(
          
          model_state ==
            "wrong",
          
          estimator ==
            "OLS",
          
          diagnostic !=
            "none"
          
        ) |>
        
        group_by(
          across(
            all_of(
              detection_group_vars
            )
          )
        ) |>
        
        summarise(
          
          null_cutoff =
            stats::quantile(
              standardized_delta,
              probs = 0.95,
              na.rm = TRUE,
              names = FALSE
            ),
          
          .groups =
            "drop"
        )
      
      
      null_parts[[j]] <-
        z
    }
    
    
    rm(x)
    
    
    if (
      j %% 250L ==
      0L
    ) {
      
      gc(
        verbose = FALSE
      )
    }
  }
  
  
  null_table <-
    dplyr::bind_rows(
      null_parts
    )
  
  
  rm(
    null_parts
  )
  
  
  gc(
    verbose = FALSE
  )
  
  
  # ----------------------------------------------------------------------------
  # PASS 2:
  # Calculate detection power for every severity level
  # ----------------------------------------------------------------------------
  
  power_parts <-
    vector(
      "list",
      length(checkpoint_files)
    )
  
  
  for (
    j in seq_along(checkpoint_files)
  ) {
    
    x <-
      readRDS(
        checkpoint_files[j]
      )
    
    
    z <-
      x |>
      
      filter(
        
        model_state ==
          "wrong",
        
        estimator ==
          "OLS",
        
        diagnostic !=
          "none"
        
      ) |>
      
      left_join(
        null_table,
        by =
          detection_group_vars
      ) |>
      
      group_by(
        
        across(
          all_of(
            c(
              detection_group_vars,
              "severity_level",
              "severity_target"
            )
          )
        )
        
      ) |>
      
      summarise(
        
        detection_power =
          mean(
            standardized_delta >
              null_cutoff,
            na.rm = TRUE
          ),
        
        .groups =
          "drop"
      )
    
    
    power_parts[[j]] <-
      z
    
    
    rm(
      x,
      z
    )
    
    
    if (
      j %% 250L ==
      0L
    ) {
      
      gc(
        verbose = FALSE
      )
    }
  }
  
  
  power_table <-
    dplyr::bind_rows(
      power_parts
    )
  
  
  rm(
    power_parts,
    null_table
  )
  
  
  gc(
    verbose = FALSE
  )
  
  
  safe_save_rds(
    power_table,
    power_path
  )
  
  
  message(
    "Detection-power table saved: ",
    power_path
  )
  
  
  # ============================================================================
  # 14. Detection boundary
  #
  # power_table is small enough to use the existing project helper normally.
  # ============================================================================
  
  
  boundary_table <-
    estimate_detection_boundary05(
      
      power_table,
      
      target_power =
        0.80
    )
  
  
  safe_save_rds(
    boundary_table,
    boundary_path
  )
  
  
  message(
    "Detection-boundary table saved: ",
    boundary_path
  )
  
  # ============================================================================
  # 15. Final validation messages
  # ============================================================================
  
  
  message(
    "Raw checkpoint files: ",
    length(checkpoint_files)
  )
  
  
  message(
    "Raw formal result rows across checkpoints: ",
    format(
      sum(
        checkpoint_row_counts
      ),
      big.mark = ",",
      scientific = FALSE
    )
  )
  
  
  message(
    "Summary rows: ",
    nrow(summary_table)
  )
  
  
  message(
    "Detection-power rows: ",
    nrow(power_table)
  )
  
  
  message(
    "Detection-boundary rows: ",
    nrow(boundary_table)
  )
  
  
  message(
    "Summary file exists: ",
    file.exists(summary_path)
  )
  
  
  message(
    "Power file exists: ",
    file.exists(power_path)
  )
  
  
  message(
    "Boundary file exists: ",
    file.exists(boundary_path)
  )
  
  
  message(
    paste0(
      "Raw formal results remain in the ",
      length(checkpoint_files),
      " checkpoint RDS files."
    )
  )
  
  
  message(
    "Script 05 completed successfully."
  )
  
  # ============================================================================
  # 16. Clean shutdown
  # ============================================================================
  
}, finally = {
  
  future::plan(
    future::sequential
  )
  
  
  try(
    parallel::stopCluster(
      cl
    ),
    silent = TRUE
  )
})