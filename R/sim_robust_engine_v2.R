# ==============================================================================
# File: /R/sim_robust_engine_v2.R
# Purpose: Enriched iteration engine for the robust comparison simulation.
#          Runs ten estimators, records complete method runtimes, detection
#          overlap, selected-set sizes, process diagnostics, coefficients,
#          standard errors, absolute bias, and 95% confidence-interval coverage.
#
# Estimators (10 total):
#   1. full       — OLS on contaminated data
#   2. cd         — OLS after Cook's-distance removal
#   3. lev        — OLS after leverage removal
#   4. dfb        — OLS after target-DFBETAS removal
#   5. mis_alpha  — exact Dinkelbach MIS with alpha_k
#   6. mis_oracle — exact Dinkelbach MIS with oracle k
#   7. mis_peel   — sigma-guided iterative peel v2
#   8. mis_sap    — selection-adjusted permutation multiscale MIS
#   9. mm         — MM-estimator
#  10. lts        — least trimmed squares
#
# Dependencies:
#   dgp_factory.R, influence_injector.R, diagnostics_classical.R,
#   estimators_robust.R, dynamic_k_adaptive.R, dinkelbach_topk.R,
#   leverage_k.R, iterative_peel_v2.R, iterative_peel_sap.R,
#   helpers_local.R
# ==============================================================================


#' Check 95% Wald Interval Coverage
#'
#' @param coef Numeric; point estimate.
#' @param se Numeric; standard error.
#' @param true_b Numeric; true coefficient value.
#'
#' @return Integer: 1L if covered, 0L if not covered, or NA_integer_.
#' @keywords internal
check_coverage_v2 <- function(coef, se, true_b) {
  coef <- unname(coef)
  se <- unname(se)
  
  if (length(coef) != 1L || length(se) != 1L ||
      !is.finite(coef) || !is.finite(se)) {
    return(NA_integer_)
  }
  
  lo <- coef - 1.96 * se
  hi <- coef + 1.96 * se
  
  as.integer(true_b >= lo && true_b <= hi)
}


#' Compute Detection Overlap Against Injected Outliers
#'
#' @param detected Integer vector of detected indices.
#' @param true_idx Integer vector of injected indices, or NULL for clean data.
#'
#' @return Fraction of injected observations recovered, or NA for clean data.
#' @keywords internal
compute_overlap <- function(detected, true_idx) {
  if (is.null(true_idx) || length(true_idx) == 0L ||
      (length(true_idx) == 1L && is.na(true_idx[1L]))) {
    return(NA_real_)
  }
  
  detected <- detected[is.finite(detected)]
  
  if (length(detected) == 0L) {
    return(0)
  }
  
  length(intersect(detected, true_idx)) / length(true_idx)
}


#' Run a Single-Shot MIS in Both Directions
#'
#' Exact Dinkelbach selection is run in both coefficient directions. The
#' direction whose cleaned OLS coefficient is closest to the MM anchor is kept.
#'
#' @param mod_full Fitted contaminated-data OLS model.
#' @param formula Model formula.
#' @param data Full data.frame.
#' @param k_val Integer; selected set size.
#' @param target_pos Integer; target coefficient position.
#' @param beta_anchor Numeric; MM reference coefficient.
#' @param res_full Named vector c(coef, se); fallback result.
#'
#' @return List with cleaned result, selected indices, and direction.
#' @keywords internal
run_mis_directional <- function(mod_full, formula, data, k_val,
                                target_pos, beta_anchor, res_full) {
  if (!is.finite(k_val) || k_val <= 0L) {
    return(list(
      result = res_full,
      indices = integer(0),
      direction = 0L
    ))
  }
  
  k_val <- as.integer(k_val)
  
  idx_pos <- dinkelbach_topk_lm(
    mod = mod_full,
    pos = target_pos,
    sign = 1L,
    k = k_val
  )
  
  idx_neg <- dinkelbach_topk_lm(
    mod = mod_full,
    pos = target_pos,
    sign = -1L,
    k = k_val
  )
  
  r_pos <- fit_clean_ols(
    formula,
    data = data,
    exclude_idx = idx_pos
  )
  
  r_neg <- fit_clean_ols(
    formula,
    data = data,
    exclude_idx = idx_neg
  )
  
  d_pos <- abs(unname(r_pos["coef"]) - beta_anchor)
  d_neg <- abs(unname(r_neg["coef"]) - beta_anchor)
  
  if (!is.finite(d_pos) && !is.finite(d_neg)) {
    return(list(
      result = res_full,
      indices = integer(0),
      direction = 0L
    ))
  }
  
  if (!is.finite(d_pos)) {
    return(list(
      result = r_neg,
      indices = idx_neg,
      direction = -1L
    ))
  }
  
  if (!is.finite(d_neg)) {
    return(list(
      result = r_pos,
      indices = idx_pos,
      direction = 1L
    ))
  }
  
  if (d_pos <= d_neg) {
    list(
      result = r_pos,
      indices = idx_pos,
      direction = 1L
    )
  } else {
    list(
      result = r_neg,
      indices = idx_neg,
      direction = -1L
    )
  }
}


#' Single Enriched Iteration of the Robust Comparison Simulation
#'
#' @param iter Integer; Monte Carlo iteration index.
#' @param n Integer; sample size.
#' @param p Integer; number of predictors.
#' @param x_type Character; design distribution.
#' @param error_type Character; error distribution.
#' @param outlier_method Character; contamination topology.
#' @param k Integer; number of injected outliers.
#' @param magnitude Numeric; contamination magnitude.
#' @param sap_alpha Numeric; SAP global-test significance level.
#' @param sap_B_perm Integer; number of residual permutations.
#' @param sap_k_grid Integer vector of candidate MIS sizes.
#' @param sap_max_iter Integer; maximum accepted SAP peel steps.
#'
#' @return A flat one-row data.frame.
#' @export
run_robust_comparison_iter_v2 <- function(
    iter,
    n = 1000L,
    p = 1L,
    x_type = "normal",
    error_type = "normal",
    outlier_method,
    k,
    magnitude,
    sap_alpha = 0.05,
    sap_B_perm = 199L,
    sap_k_grid = c(1L, 2L, 5L, 10L, 20L, 50L, 100L),
    sap_max_iter = 1L
) {
  
  # =================================================================
  # 1. Data generation and contamination
  # =================================================================
  dat_clean <- generate_complex_data(
    n = n,
    p = p,
    x_type = x_type,
    error_type = error_type
  )
  
  true_b <- dat_clean$true_beta[1L]
  
  if (outlier_method != "none") {
    dat <- apply_influence_shift(
      dat_clean,
      method = outlier_method,
      k = k,
      magnitude = magnitude
    )
    true_idx <- dat$outlier_indices
  } else {
    dat <- dat_clean
    true_idx <- NULL
  }
  
  df <- data.frame(
    y = dat$y,
    x = dat$X[, 1L]
  )
  
  # =================================================================
  # 2. Full OLS
  # =================================================================
  t0 <- proc.time()[3L]
  
  mod_full <- stats::lm(y ~ x, data = df)
  res_full <- fit_clean_ols(
    y ~ x,
    data = df,
    exclude_idx = integer(0)
  )
  
  cpu_full <- proc.time()[3L] - t0
  
  # =================================================================
  # 3. Shared MM anchor for alpha_k and direction selection
  # =================================================================
  t0 <- proc.time()[3L]
  
  mod_mm_obj <- tryCatch(
    robustbase::lmrob(
      y ~ x,
      data = df,
      setting = "KS2014"
    ),
    error = function(e) NULL
  )
  
  cpu_mm_anchor <- proc.time()[3L] - t0

  # Validate the shared MM fit.
  mm_scale <- if (
    !is.null(mod_mm_obj) &&
    length(mod_mm_obj$scale) == 1L
  ) {
    unname(mod_mm_obj$scale)
  } else {
    NA_real_
  }
  
  mm_slope <- if (!is.null(mod_mm_obj)) {
    tryCatch(
      unname(stats::coef(mod_mm_obj)["x"]),
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }
  
  mm_valid <- (
    !is.null(mod_mm_obj) &&
      isTRUE(mod_mm_obj$converged) &&
      is.finite(mm_scale) &&
      mm_scale > sqrt(.Machine$double.eps) &&
      is.finite(mm_slope)
  )
  
  mm_zero_scale <- (
    !is.null(mod_mm_obj) &&
      is.finite(mm_scale) &&
      mm_scale <= sqrt(.Machine$double.eps)
  )
  
  # Use MM only when the fit is numerically valid.
  beta_mm <- if (mm_valid) {
    mm_slope
  } else {
    unname(stats::coef(mod_full)["x"])
  }
  
  # Preserve the existing output-column name used by Script 04.
  mm_converged <- mm_valid
  
  # =================================================================
  # 4. Classical diagnostics
  # Complete runtime: detection + cleaned OLS refit
  # =================================================================
  t0 <- proc.time()[3L]
  
  cd_idx <- get_classical_set(
    mod_full,
    target_var = "x",
    k = NULL,
    metric = "cooks_d"
  )
  res_cd <- fit_clean_ols(
    y ~ x,
    data = df,
    exclude_idx = cd_idx
  )
  
  cpu_cd <- proc.time()[3L] - t0
  
  t0 <- proc.time()[3L]
  
  lev_idx <- get_classical_set(
    mod_full,
    target_var = "x",
    k = NULL,
    metric = "leverage"
  )
  res_lev <- fit_clean_ols(
    y ~ x,
    data = df,
    exclude_idx = lev_idx
  )
  
  cpu_lev <- proc.time()[3L] - t0
  
  t0 <- proc.time()[3L]
  
  dfb_idx <- get_classical_set(
    mod_full,
    target_var = "x",
    k = NULL,
    metric = "dfbetas_target"
  )
  res_dfb <- fit_clean_ols(
    y ~ x,
    data = df,
    exclude_idx = dfb_idx
  )
  
  cpu_dfb <- proc.time()[3L] - t0
  
  # =================================================================
  # 5. MIS-alpha
  # Complete runtime: MM anchor + k selection + search + refit
  # =================================================================
  t0 <- proc.time()[3L]
  
  k_alpha_val <- if (mm_valid) {
    tryCatch(
      as.integer(alpha_k(mod_mm_obj)),
      error = function(e) 0L
    )
  } else {
    0L
  }
  
  mis_alpha <- run_mis_directional(
    mod_full = mod_full,
    formula = y ~ x,
    data = df,
    k_val = k_alpha_val,
    target_pos = 2L,
    beta_anchor = beta_mm,
    res_full = res_full
  )
  
  cpu_mis_alpha <- cpu_mm_anchor + (proc.time()[3L] - t0)
  
  # =================================================================
  # 6. MIS-oracle
  # Complete runtime: MM anchor + oracle k + search + refit
  # =================================================================
  t0 <- proc.time()[3L]
  
  k_oracle_val <- oracle_k(
    if (outlier_method == "none") 0L else k
  )
  
  mis_oracle <- run_mis_directional(
    mod_full = mod_full,
    formula = y ~ x,
    data = df,
    k_val = k_oracle_val,
    target_pos = 2L,
    beta_anchor = beta_mm,
    res_full = res_full
  )
  
  cpu_mis_oracle <- cpu_mm_anchor + (proc.time()[3L] - t0)
  
  # =================================================================
  # 7. Sigma-guided iterative peel v2
  # Complete runtime: peeling + final cleaned OLS refit
  # =================================================================
  t0 <- proc.time()[3L]
  
  peel_v2_result <- tryCatch(
    iterative_peel_v2(
      formula = y ~ x,
      data = df,
      target_var = "x",
      target_pos = 2L,
      batch_size = 1L,
      max_iter = 50L,
      max_k_frac = 0.06,
      detector = "dinkelbach",
      k_method = "leverage",
      verbose = FALSE
    ),
    error = function(e) {
      list(
        excluded = integer(0),
        k_total = 0L,
        n_iters = 0L,
        stop_reason = "error",
        beta_trajectory = numeric(0),
        sigma_trajectory = numeric(0),
        error_message = conditionMessage(e)
      )
    }
  )
  
  res_peel_v2 <- fit_clean_ols(
    y ~ x,
    data = df,
    exclude_idx = peel_v2_result$excluded
  )
  
  cpu_peel_v2 <- proc.time()[3L] - t0
  
  # =================================================================
  # 8. Selection-Adjusted Permutation MIS
  # Complete runtime: permutation search + final cleaned OLS refit
  # =================================================================
  t0 <- proc.time()[3L]
  
  peel_sap_result <- tryCatch(
    iterative_peel_sap(
      formula = y ~ x,
      data = df,
      target_var = "x",
      target_pos = 2L,
      k_grid = sap_k_grid,
      B_perm = sap_B_perm,
      alpha = sap_alpha,
      max_iter = sap_max_iter,
      max_k_frac = 0.10,
      verbose = FALSE
    ),
    error = function(e) {
      list(
        excluded = integer(0),
        k_total = 0L,
        n_iters = 0L,
        stop_reason = "error",
        beta_trajectory = numeric(0),
        sigma_trajectory = numeric(0),
        global_p_trajectory = NA_real_,
        selected_k_trajectory = 0L,
        direction_trajectory = 0L,
        excess_ratio_trajectory = NA_real_,
        profile_trajectory = list(),
        error_message = conditionMessage(e)
      )
    }
  )
  
  res_peel_sap <- fit_clean_ols(
    y ~ x,
    data = df,
    exclude_idx = peel_sap_result$excluded
  )
  
  cpu_peel_sap <- proc.time()[3L] - t0
  
  sap_p_values <- peel_sap_result$global_p_trajectory
  sap_p_values <- sap_p_values[is.finite(sap_p_values)]
  
  sap_final_p <- if (length(sap_p_values) > 0L) {
    tail(sap_p_values, 1L)
  } else {
    NA_real_
  }
  
  sap_min_p <- if (length(sap_p_values) > 0L) {
    min(sap_p_values)
  } else {
    NA_real_
  }
  
  sap_directions <- peel_sap_result$direction_trajectory
  sap_directions <- sap_directions[is.finite(sap_directions)]
  
  sap_direction <- if (length(sap_directions) > 0L) {
    as.integer(tail(sap_directions, 1L))
  } else {
    0L
  }
  
  sap_excess <- peel_sap_result$excess_ratio_trajectory
  sap_excess <- sap_excess[is.finite(sap_excess)]
  
  sap_peak_excess <- if (length(sap_excess) > 0L) {
    max(sap_excess)
  } else {
    NA_real_
  }
  
  # =================================================================
  # 9. Direct robust estimators
  # =================================================================
  t0 <- proc.time()[3L]
  
  res_mm <- fit_mm_estimator(
    y ~ x,
    data = df
  )
  
  cpu_mm <- proc.time()[3L] - t0
  
  t0 <- proc.time()[3L]
  
  res_lts <- fit_lts_estimator(
    y ~ x,
    data = df
  )
  
  cpu_lts <- proc.time()[3L] - t0
  
  # =================================================================
  # 10. Detection overlap
  # =================================================================
  overlap_cd <- compute_overlap(cd_idx, true_idx)
  overlap_lev <- compute_overlap(lev_idx, true_idx)
  overlap_dfb <- compute_overlap(dfb_idx, true_idx)
  overlap_mis_alpha <- compute_overlap(mis_alpha$indices, true_idx)
  overlap_mis_oracle <- compute_overlap(mis_oracle$indices, true_idx)
  overlap_peel_v2 <- compute_overlap(peel_v2_result$excluded, true_idx)
  overlap_peel_sap <- compute_overlap(peel_sap_result$excluded, true_idx)
  
  # =================================================================
  # 11. Assemble flat output
  # =================================================================
  data.frame(
    iter = iter,
    x_type = x_type,
    error_type = error_type,
    outlier_method = outlier_method,
    set_size = if (outlier_method == "none") 0L else k,
    
    k_cd = length(cd_idx),
    k_lev = length(lev_idx),
    k_dfb = length(dfb_idx),
    k_alpha = k_alpha_val,
    k_oracle = k_oracle_val,
    k_peel_v2 = peel_v2_result$k_total,
    k_peel_sap = peel_sap_result$k_total,
    
    overlap_cd = overlap_cd,
    overlap_lev = overlap_lev,
    overlap_dfb = overlap_dfb,
    overlap_mis_alpha = overlap_mis_alpha,
    overlap_mis_oracle = overlap_mis_oracle,
    overlap_peel_v2 = overlap_peel_v2,
    overlap_peel_sap = overlap_peel_sap,
    
    dir_alpha = mis_alpha$direction,
    dir_oracle = mis_oracle$direction,
    
    peel_v2_stop = peel_v2_result$stop_reason,
    peel_v2_iters = peel_v2_result$n_iters,
    
    peel_sap_stop = peel_sap_result$stop_reason,
    peel_sap_iters = peel_sap_result$n_iters,
    peel_sap_final_p = sap_final_p,
    peel_sap_min_p = sap_min_p,
    peel_sap_direction = sap_direction,
    peel_sap_peak_excess = sap_peak_excess,
    
    mm_converged = mm_converged,
    mm_valid = mm_valid,
    mm_zero_scale = mm_zero_scale,
    mm_scale = mm_scale,
    
    coef_full = unname(res_full["coef"]),
    coef_cd = unname(res_cd["coef"]),
    coef_lev = unname(res_lev["coef"]),
    coef_dfb = unname(res_dfb["coef"]),
    coef_mis_alpha = unname(mis_alpha$result["coef"]),
    coef_mis_oracle = unname(mis_oracle$result["coef"]),
    coef_mis_peel = unname(res_peel_v2["coef"]),
    coef_mis_sap = unname(res_peel_sap["coef"]),
    coef_mm = unname(res_mm["coef"]),
    coef_lts = unname(res_lts["coef"]),
    
    se_full = unname(res_full["se"]),
    se_cd = unname(res_cd["se"]),
    se_lev = unname(res_lev["se"]),
    se_dfb = unname(res_dfb["se"]),
    se_mis_alpha = unname(mis_alpha$result["se"]),
    se_mis_oracle = unname(mis_oracle$result["se"]),
    se_mis_peel = unname(res_peel_v2["se"]),
    se_mis_sap = unname(res_peel_sap["se"]),
    se_mm = unname(res_mm["se"]),
    se_lts = unname(res_lts["se"]),
    
    bias_full = unname(abs(res_full["coef"] - true_b)),
    bias_cd = unname(abs(res_cd["coef"] - true_b)),
    bias_lev = unname(abs(res_lev["coef"] - true_b)),
    bias_dfb = unname(abs(res_dfb["coef"] - true_b)),
    bias_mis_alpha = unname(abs(mis_alpha$result["coef"] - true_b)),
    bias_mis_oracle = unname(abs(mis_oracle$result["coef"] - true_b)),
    bias_mis_peel = unname(abs(res_peel_v2["coef"] - true_b)),
    bias_mis_sap = unname(abs(res_peel_sap["coef"] - true_b)),
    bias_mm = unname(abs(res_mm["coef"] - true_b)),
    bias_lts = unname(abs(res_lts["coef"] - true_b)),
    
    cov_full = check_coverage_v2(
      res_full["coef"], res_full["se"], true_b
    ),
    cov_cd = check_coverage_v2(
      res_cd["coef"], res_cd["se"], true_b
    ),
    cov_lev = check_coverage_v2(
      res_lev["coef"], res_lev["se"], true_b
    ),
    cov_dfb = check_coverage_v2(
      res_dfb["coef"], res_dfb["se"], true_b
    ),
    cov_mis_alpha = check_coverage_v2(
      mis_alpha$result["coef"], mis_alpha$result["se"], true_b
    ),
    cov_mis_oracle = check_coverage_v2(
      mis_oracle$result["coef"], mis_oracle$result["se"], true_b
    ),
    cov_mis_peel = check_coverage_v2(
      res_peel_v2["coef"], res_peel_v2["se"], true_b
    ),
    cov_mis_sap = check_coverage_v2(
      res_peel_sap["coef"], res_peel_sap["se"], true_b
    ),
    cov_mm = check_coverage_v2(
      res_mm["coef"], res_mm["se"], true_b
    ),
    cov_lts = check_coverage_v2(
      res_lts["coef"], res_lts["se"], true_b
    ),
    
    cpu_full = unname(cpu_full),
    cpu_cd = unname(cpu_cd),
    cpu_lev = unname(cpu_lev),
    cpu_dfb = unname(cpu_dfb),
    cpu_mis_alpha = unname(cpu_mis_alpha),
    cpu_mis_oracle = unname(cpu_mis_oracle),
    cpu_peel_v2 = unname(cpu_peel_v2),
    cpu_peel_sap = unname(cpu_peel_sap),
    cpu_mm = unname(cpu_mm),
    cpu_lts = unname(cpu_lts),
    
    stringsAsFactors = FALSE
  )
}
