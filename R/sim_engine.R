# ==============================================================================
# File: /R/sim_engine.R
# Purpose: Core simulation wrapper combining DGP, influence injection, MIS 
#          detection, classical LOO diagnostics, and EVD parameter estimation. 
#          Orchestrates one complete Monte Carlo iteration: generates data under 
#          a specified DGP, optionally injects adversarial influence, detects the 
#          Most Influential Set via exact sensitivity search, computes classical 
#          diagnostics (Cook's D, Leverage, DFBETAS) as baselines, and fits the 
#          block-maxima Extreme Value Distribution.
# Dependencies: Requires /R/dgp_factory.R, /R/influence_injector.R, 
#               /R/evt_iter_dm.R, and /R/exact_dfb_bmx.R to be sourced.
# ==============================================================================
#'
#' @param iter Integer: Iteration ID for tracking across simulation grid
#' @param n Integer: Sample size for the generated dataset
#' @param p Integer: Number of predictors (currently only p = 1 supported)
#' @param x_type Character: Distribution type for the design matrix X 
#'        (from dgp_factory: "normal", "skewed_t", "pareto")
#' @param error_type Character: Distribution type for the error term 
#'        (from dgp_factory: "normal", "mixed_normal", "skewed_t", "golm",
#'        "beta_logistic", "gpd", "contaminated", "pareto")
#' @param outlier_method Character: Injection strategy from influence_injector
#'        ("none", "vertical_outlier", "good_leverage", "bad_leverage")
#' @param k Integer: Size of the influential set to inject and detect (default = 1)
#' @param magnitude Numeric: Scale multiplier for the adversarial shift (must be > 0, default = 5)
#' @param block_count Integer: Number of blocks for block-maxima EVD estimation (default = 40)
#' @param dist_param Numeric: Shape/tail parameter passed to dgp_factory (default = 3)
#' @param mix_prop Numeric: Mixture proportion passed to dgp_factory (default = 0.1)
#'
#' @return A flat, 1-row data.frame containing all simulation inputs, detection
#'         metrics, EVD parameters (shape, scale, loc), the observed set DFBETA,
#'         the extreme value p-value, convergence flag, and wall-clock time.
#' @export
run_mis_iteration <- function(iter = 1, n = 1000, p = 1,
                              x_type = "normal", error_type = "normal",
                              outlier_method = "none", k = 1, magnitude = 5,
                              block_count = 40, dist_param = 3, mix_prop = 0.1) {
  
  if (p != 1) stop("sim_engine currently only supports p = 1.")
  start_time <- Sys.time()
  
  # -------------------------------------------------------------------------
  # 1. Data Generation & Injection
  # -------------------------------------------------------------------------
  dat_clean <- generate_complex_data(n = n, p = p, x_type = x_type,
                                     error_type = error_type,
                                     dist_param = dist_param,
                                     mix_prop = mix_prop)
  
  if (outlier_method != "none") {
    dat <- apply_influence_shift(dat_clean, method = outlier_method,
                                 k = k, magnitude = magnitude)
    true_injected_indices <- dat$outlier_indices
  } else {
    dat <- dat_clean
    true_injected_indices <- NA
  }
  
  df_sim <- data.frame(y = dat$y, x = dat$X[, 1])
  
  # -------------------------------------------------------------------------
  # 2. Fit Base Model & Detect MIS
  # -------------------------------------------------------------------------
  base_model <- lm(y ~ x, data = df_sim)
  
  target_sign <- sign(coef(base_model)[2])
  if (is.na(target_sign) || target_sign == 0) target_sign <- 1
  
  sens_obj_pos <- suppressWarnings({
    influence::sens(base_model,
                    lambda = influence::set_lambda("beta_i", pos = 2, sign = 1))
  })
  sens_obj_neg <- suppressWarnings({
    influence::sens(base_model,
                    lambda = influence::set_lambda("beta_i", pos = 2, sign = -1))
  })
  
  empirical_mis_pos <- sens_obj_pos$influence$id[1:k]
  empirical_mis_neg <- sens_obj_neg$influence$id[1:k]
  
  if (outlier_method != "none") {
    # Pick the MIS direction with better overlap for detection metrics
    overlap_pos <- length(intersect(true_injected_indices, empirical_mis_pos))
    overlap_neg <- length(intersect(true_injected_indices, empirical_mis_neg))
    if (overlap_neg > overlap_pos) {
      empirical_mis <- empirical_mis_neg
    } else {
      empirical_mis <- empirical_mis_pos
    }
    detection_success <- (length(intersect(true_injected_indices, 
                                           empirical_mis)) / k) >= 0.80
  } else {
    detection_success <- NA
    empirical_mis <- empirical_mis_pos
  }
  
  # 2b. Classical LOO Diagnostics
  
  cooks_vals   <- cooks.distance(base_model)
  lev_vals     <- hatvalues(base_model)
  dfbetas_vals <- dfbetas(base_model)[, "x"]  # slope coefficient
  
  top_k_cooks   <- order(abs(cooks_vals),   decreasing = TRUE)[1:k]
  top_k_lev     <- order(abs(lev_vals),     decreasing = TRUE)[1:k]
  top_k_dfbetas <- order(abs(dfbetas_vals), decreasing = TRUE)[1:k]
  
  if (outlier_method != "none") {
    detect_cooks   <- (length(intersect(true_injected_indices, 
                                        top_k_cooks)) / k) >= 0.80
    detect_lev     <- (length(intersect(true_injected_indices, 
                                        top_k_lev)) / k) >= 0.80
    detect_dfbetas <- (length(intersect(true_injected_indices, 
                                        top_k_dfbetas)) / k) >= 0.80
    
    overlap_mis     <- length(intersect(true_injected_indices, empirical_mis)) / k
    overlap_cooks   <- length(intersect(true_injected_indices, top_k_cooks)) / k
    overlap_lev     <- length(intersect(true_injected_indices, top_k_lev)) / k
    overlap_dfbetas <- length(intersect(true_injected_indices, top_k_dfbetas)) / k
  } else {
    detect_cooks   <- NA
    detect_lev     <- NA
    detect_dfbetas <- NA
    overlap_mis     <- NA
    overlap_cooks   <- NA
    overlap_lev     <- NA
    overlap_dfbetas <- NA
  }
  
  # -------------------------------------------------------------------------
  # 3. Estimate Extreme Value Distribution (EVD) Parameters
  # -------------------------------------------------------------------------
  # POST-SELECTION FIX: Choose the correct set for the EVD test.
  #   The EVD null is calibrated for a PRE-SPECIFIED (fixed) candidate set.
  #   Testing the MIS (selected to maximise influence) causes 100% rejection
  #   even on clean data — classic post-selection bias.
  #
  #   - With outliers: test the KNOWN INJECTED indices (ground truth).
  #     This answers: "is the injected set's influence excessive under the EVD?"
  #   - Without outliers: test a RANDOM set of size k.
  #     This answers: "does the EVD test maintain nominal size?"
  
  if (outlier_method != "none") {
    evd_test_set <- true_injected_indices
  } else {
    evd_test_set <- sample(seq_len(n), k)
  }
  
  Z_matrix <- matrix(1, nrow = length(dat$y), ncol = 1)
  
  res_exact <- tryCatch({
    evt_iter_dm(
      y = dat$y,
      x = dat$X[, 1],
      Z = Z_matrix,
      set = evd_test_set,
      block_count = block_count
    )
  }, error = function(e) {
    warning(sprintf("evt_iter_dm failed (iter=%d, x=%s, err=%s, outlier=%s): %s",
                    iter, x_type, error_type, outlier_method,
                    conditionMessage(e)))
    data.frame(shape = NA_real_, scale = NA_real_, loc = NA_real_,
               set_dfb = NA_real_, p_value = NA_real_, converged = FALSE,
               stringsAsFactors = FALSE)
  })
  
  # If evt_iter_dm returned a list, coerce to 1-row data.frame
  expected_names <- c("shape", "scale", "loc", "set_dfb", "p_value", "converged")
  if (is.list(res_exact) && !is.data.frame(res_exact)) {
    missing <- setdiff(expected_names, names(res_exact))
    if (length(missing) > 0) {
      stop(sprintf("evt_iter_dm list is missing fields: %s  (got: %s)",
                   paste(missing, collapse = ", "),
                   paste(names(res_exact), collapse = ", ")))
    }
    res_exact <- as.data.frame(res_exact[expected_names],
                               stringsAsFactors = FALSE)
  }
  
  end_time <- Sys.time()
  compute_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # -------------------------------------------------------------------------
  # 4. Assemble Tidy Output
  # -------------------------------------------------------------------------
  res_df <- data.frame(
    iter = iter,
    n_obs = n,
    x_type = x_type,
    error_type = error_type,
    dist_param = dist_param,
    outlier_method = outlier_method,
    magnitude = magnitude,
    set_size = k,
    block_count = block_count,
    # MIS detection
    detection_success = detection_success,
    # Classical detection (exact set match)
    detect_cooks = detect_cooks,
    detect_lev = detect_lev,
    detect_dfbetas = detect_dfbetas,
    # Partial overlap (fraction of injected points found in top-k)
    overlap_mis = overlap_mis,
    overlap_cooks = overlap_cooks,
    overlap_lev = overlap_lev,
    overlap_dfbetas = overlap_dfbetas,
    # Which set was tested by EVD
    evd_test_type = if (outlier_method != "none") "injected" else "random",
    # Timing
    compute_time = compute_time,
    stringsAsFactors = FALSE
  )
  
  if (!is.data.frame(res_exact) || nrow(res_exact) != 1) {
    stop("evt_iter_dm returned unexpected output structure.")
  }
  
  res_final <- cbind(res_df, res_exact)
  return(res_final)
}