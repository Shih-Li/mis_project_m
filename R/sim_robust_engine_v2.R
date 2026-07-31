# ==============================================================================
# File: /R/sim_robust_engine_v2.R
# Purpose: Enriched iteration engine for the robust comparison simulation.
#          Runs eleven estimators/tools, records complete method runtimes, 
#          overlap, selected-set sizes, process diagnostics, coefficients,
#          standard errors, absolute bias, and 95% confidence-interval coverage.
#
# Estimators/tools (11 total):
#   1. full       — OLS on contaminated data
#   2. cd         — OLS after Cook's-distance removal
#   3. lev        — OLS after leverage removal
#   4. dfb        — OLS after target-DFBETAS removal
#   5. mis_alpha  — Dinkelbach MIS with alpha_k
#   6. mis_oracle — Dinkelbach MIS with oracle k
#   7. mis_peel   — sigma-guided iterative peel v2
#   8. mis_sap    — formal selection-adjusted permutation MIS
#   9. mis_sek    — independent effective-k certification tool
#  10. mm         — MM-estimator
#  11. lts        — least trimmed squares
#
# Dependencies:
#   dgp_factory.R, influence_injector.R, diagnostics_classical.R,
#   estimators_robust.R, dynamic_k_adaptive.R, dinkelbach_topk.R,
#   leverage_k.R, iterative_peel_v2.R, mis_sap.R, mis_sek.R,
#   mis_sek_adapter.R, helpers_local.R
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

#' Collapse an integer set for flat simulation output
#' @keywords internal
collapse_integer_set_v2 <- function(x) {
  if (is.null(x)) {
    return(NA_character_)
  }
  
  if (length(x) == 0L) {
    return("")
  }
  
  paste(
    sort(unique(as.integer(x))),
    collapse = "|"
  )
}


#' Collapse character diagnostics for flat simulation output
#' @keywords internal
collapse_character_set_v2 <- function(x) {
  if (is.null(x)) {
    return(NA_character_)
  }
  
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  
  if (length(x) == 0L) {
    return("")
  }
  
  paste(unique(x), collapse = "|")
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
#' @param sap_max_fraction Numeric; maximum SAP coalition fraction.
#' @param sek_alpha Numeric; MIS-sek calibration level.
#' @param sek_B_cal Integer; independent MIS-sek calibration repetitions.
#' @param sek_k_grid Integer vector of MIS-sek candidate sizes.
#' @param sek_max_fraction Numeric; maximum MIS-sek candidate fraction.
#' @param sek_eta_n Numeric; MIS-sek near-optimality slack.
#' @param sek_minimum_denominator_fraction Numeric; denominator safety floor.
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
    sap_max_fraction = 0.05,
    
    sek_alpha = 0.05,
    sek_B_cal = 199L,
    sek_k_grid = c(1L, 2L, 5L, 10L, 20L, 50L, 100L),
    sek_max_fraction = 0.05,
    sek_eta_n = 0,
    sek_minimum_denominator_fraction = 0.05
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
  # 8. Formal MIS-SAP
  # Complete runtime: calibrated search + optional cleaning + OLS refit
  # =================================================================
  t0 <- proc.time()[3L]
  
  sap_run <- tryCatch(
    list(
      result = mis_sap(
        formula = y ~ x,
        data = df,
        target = "x",
        k_grid = sap_k_grid,
        B_perm = sap_B_perm,
        alpha = sap_alpha,
        max_fraction = sap_max_fraction,
        use_robust_anchor = TRUE,
        anchor_coefficient = NULL
      ),
      error = NA_character_
    ),
    error = function(e) {
      list(
        result = NULL,
        error = conditionMessage(e)
      )
    }
  )
  
  mis_sap_result <- sap_run$result
  mis_sap_error <- sap_run$error
  mis_sap_ok <- inherits(
    mis_sap_result,
    "mis_sap"
  )
  
  sap_sensitivity_idx <- if (
    mis_sap_ok &&
    isTRUE(mis_sap_result$global$reject)
  ) {
    as.integer(mis_sap_result$sensitivity$set)
  } else {
    integer(0L)
  }
  
  sap_cleaning_idx <- if (
    mis_sap_ok &&
    isTRUE(mis_sap_result$cleaning$permitted)
  ) {
    as.integer(mis_sap_result$cleaning$set)
  } else {
    integer(0L)
  }
  
  # Only the formal cleaning recommendation is used for the estimator.
  # No detection or no permitted cleaning means no deletion.
  res_mis_sap <- if (mis_sap_ok) {
    fit_clean_ols(
      y ~ x,
      data = df,
      exclude_idx = sap_cleaning_idx
    )
  } else {
    c(
      coef = NA_real_,
      se = NA_real_
    )
  }
  
  sap_state <- if (mis_sap_ok) {
    as.character(mis_sap_result$state)
  } else {
    "error"
  }
  
  sap_reject <- if (mis_sap_ok) {
    isTRUE(mis_sap_result$global$reject)
  } else {
    NA
  }
  
  sap_global_p <- if (mis_sap_ok) {
    unname(mis_sap_result$global$p_value)
  } else {
    NA_real_
  }
  
  sap_stop_reason <- if (mis_sap_ok) {
    as.character(mis_sap_result$global$stop_reason)
  } else {
    "error"
  }
  
  sap_sensitivity_k <- if (mis_sap_ok) {
    as.integer(mis_sap_result$sensitivity$selected_k)
  } else {
    NA_integer_
  }
  
  sap_sensitivity_direction <- if (mis_sap_ok) {
    as.integer(mis_sap_result$sensitivity$direction)
  } else {
    NA_integer_
  }
  
  sap_cleaning_permitted <- if (mis_sap_ok) {
    isTRUE(mis_sap_result$cleaning$permitted)
  } else {
    NA
  }
  
  sap_cleaning_direction <- if (mis_sap_ok) {
    as.integer(mis_sap_result$cleaning$direction)
  } else {
    NA_integer_
  }
  
  sap_peak_excess <- if (mis_sap_ok) {
    unname(mis_sap_result$global$peak_excess_ratio)
  } else {
    NA_real_
  }
  
  cpu_mis_sap <- proc.time()[3L] - t0
  
  # =================================================================
  # 9. Independent MIS-sek certification
  # =================================================================
  t0 <- proc.time()[3L]
  
  sek_run <- tryCatch(
    list(
      bundle = run_mis_sek_from_data(
        formula = y ~ x,
        data = df,
        target = "x",
        k_grid = sek_k_grid,
        B_cal = sek_B_cal,
        alpha = sek_alpha,
        max_fraction = sek_max_fraction,
        eta_n = sek_eta_n,
        minimum_denominator_fraction =
          sek_minimum_denominator_fraction
      ),
      error = NA_character_
    ),
    error = function(e) {
      list(
        bundle = NULL,
        error = conditionMessage(e)
      )
    }
  )
  
  mis_sek_bundle <- sek_run$bundle
  mis_sek_error <- sek_run$error
  
  mis_sek_ok <- (
    !is.null(mis_sek_bundle) &&
      is.list(mis_sek_bundle) &&
      inherits(
        mis_sek_bundle$certificate,
        "mis_sek"
      )
  )
  
  mis_sek_result <- if (mis_sek_ok) {
    mis_sek_bundle$certificate
  } else {
    NULL
  }
  
  sek_point_candidate <- (
    mis_sek_ok &&
      identical(
        mis_sek_result$state,
        "stable_effective_k"
      ) &&
      identical(
        mis_sek_result$selection_type,
        "point"
      ) &&
      isTRUE(
        mis_sek_result$automatic_top_k_submission
      ) &&
      length(mis_sek_result$selected_k) == 1L &&
      is.finite(mis_sek_result$selected_k)
  )
  
  sek_point_idx <- if (
    sek_point_candidate &&
    !is.null(mis_sek_bundle$point_set)
  ) {
    sort(
      unique(
        as.integer(
          mis_sek_bundle$point_set
        )
      )
    )
  } else {
    integer(0L)
  }
  
  sek_point_available <- (
    sek_point_candidate &&
      length(sek_point_idx) ==
      as.integer(mis_sek_result$selected_k) &&
      all(
        sek_point_idx >= 1L &
          sek_point_idx <= nrow(df)
      )
  )
  
  res_mis_sek <- if (sek_point_available) {
    fit_clean_ols(
      y ~ x,
      data = df,
      exclude_idx = sek_point_idx
    )
  } else {
    c(
      coef = NA_real_,
      se = NA_real_
    )
  }
  
  sek_state <- if (mis_sek_ok) {
    as.character(mis_sek_result$state)
  } else {
    "error"
  }
  
  sek_selection_type <- if (mis_sek_ok) {
    as.character(mis_sek_result$selection_type)
  } else {
    "none"
  }
  
  sek_selected_k <- if (mis_sek_ok) {
    as.integer(mis_sek_result$selected_k)
  } else {
    NA_integer_
  }
  
  sek_selected_k_set <- if (mis_sek_ok) {
    mis_sek_result$selected_k_set
  } else {
    integer(0L)
  }
  
  sek_selected_k_range <- if (mis_sek_ok) {
    mis_sek_result$selected_k_range
  } else {
    c(
      lower = NA_integer_,
      upper = NA_integer_
    )
  }
  
  sek_selected_sgn <- if (mis_sek_ok) {
    as.integer(mis_sek_result$selected_sgn)
  } else {
    NA_integer_
  }
  
  sek_global_reject <- if (mis_sek_ok) {
    isTRUE(mis_sek_result$global$reject)
  } else {
    NA
  }
  
  sek_support_complete <- if (mis_sek_ok) {
    isTRUE(
      mis_sek_result$stability$support_complete
    )
  } else {
    NA
  }
  
  sek_direction_stable <- if (mis_sek_ok) {
    isTRUE(
      mis_sek_result$stability$direction_stable
    )
  } else {
    NA
  }
  
  sek_denominator_safe <- if (mis_sek_ok) {
    isTRUE(
      mis_sek_result$stability$denominator_safe
    )
  } else {
    NA
  }
  
  sek_formal_guarantee_eligible <- if (mis_sek_ok) {
    isTRUE(
      mis_sek_result$stability$
        formal_guarantee_eligible
    )
  } else {
    NA
  }
  
  sek_instability_reasons <- if (mis_sek_ok) {
    collapse_character_set_v2(
      mis_sek_result$stability$
        instability_reasons
    )
  } else {
    NA_character_
  }
  
  cpu_mis_sek <- proc.time()[3L] - t0
  
  # =================================================================
  # 10. Direct robust estimators
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
  # 11. Detection overlap
  # =================================================================
  overlap_cd <- compute_overlap(cd_idx, true_idx)
  overlap_lev <- compute_overlap(lev_idx, true_idx)
  overlap_dfb <- compute_overlap(dfb_idx, true_idx)
  overlap_mis_alpha <- compute_overlap(mis_alpha$indices, true_idx)
  overlap_mis_oracle <- compute_overlap(mis_oracle$indices, true_idx)
  overlap_peel_v2 <- compute_overlap(peel_v2_result$excluded, true_idx)
  overlap_mis_sap_sensitivity <- compute_overlap(
    sap_sensitivity_idx,
    true_idx
  )
  
  overlap_mis_sap_cleaning <- compute_overlap(
    sap_cleaning_idx,
    true_idx
  )
  
  overlap_mis_sek <- if (sek_point_available) {
    compute_overlap(
      sek_point_idx,
      true_idx
    )
  } else {
    NA_real_
  }
  
  # =================================================================
  # 12. Assemble flat output
  # =================================================================
  data.frame(
    iter = iter,
    
    n_obs = n,
    design_k = as.integer(k),
    contam_prop = if (
      outlier_method == "none"
    ) {
      0
    } else {
      k / n
    },
    
    x_type = x_type,
    error_type = error_type,
    outlier_method = outlier_method,
    
    set_size = if (
      outlier_method == "none"
    ) {
      0L
    } else {
      as.integer(k)
    },
    
    k_cd = length(cd_idx),
    k_lev = length(lev_idx),
    k_dfb = length(dfb_idx),
    k_alpha = k_alpha_val,
    k_oracle = k_oracle_val,
    k_peel_v2 = peel_v2_result$k_total,
    
    k_mis_sap_sensitivity =
      length(sap_sensitivity_idx),
    
    k_mis_sap_cleaning =
      length(sap_cleaning_idx),
    
    k_mis_sek = if (
      sek_point_available
    ) {
      as.integer(sek_selected_k)
    } else {
      NA_integer_
    },
    
    overlap_cd = overlap_cd,
    overlap_lev = overlap_lev,
    overlap_dfb = overlap_dfb,
    overlap_mis_alpha = overlap_mis_alpha,
    overlap_mis_oracle = overlap_mis_oracle,
    overlap_peel_v2 = overlap_peel_v2,
    
    overlap_mis_sap_sensitivity =
      overlap_mis_sap_sensitivity,
    
    overlap_mis_sap_cleaning =
      overlap_mis_sap_cleaning,
    
    overlap_mis_sek =
      overlap_mis_sek,
    
    dir_alpha = mis_alpha$direction,
    dir_oracle = mis_oracle$direction,
    
    peel_v2_stop = peel_v2_result$stop_reason,
    peel_v2_iters = peel_v2_result$n_iters,
    
    mis_sap_state = sap_state,
    mis_sap_ok = mis_sap_ok,
    mis_sap_reject = sap_reject,
    mis_sap_global_p = sap_global_p,
    mis_sap_stop_reason = sap_stop_reason,
    
    mis_sap_sensitivity_k =
      sap_sensitivity_k,
    
    mis_sap_sensitivity_direction =
      sap_sensitivity_direction,
    
    mis_sap_cleaning_permitted =
      sap_cleaning_permitted,
    
    mis_sap_cleaning_direction =
      sap_cleaning_direction,
    
    mis_sap_peak_excess =
      sap_peak_excess,
    
    mis_sap_error =
      mis_sap_error,
    
    mis_sek_state =
      sek_state,
    
    mis_sek_ok =
      mis_sek_ok,
    
    mis_sek_global_reject =
      sek_global_reject,
    
    mis_sek_selection_type =
      sek_selection_type,
    
    mis_sek_point_available =
      sek_point_available,
    
    mis_sek_selected_k =
      sek_selected_k,
    
    mis_sek_selected_k_set =
      collapse_integer_set_v2(
        sek_selected_k_set
      ),
    
    mis_sek_k_lower =
      unname(
        sek_selected_k_range["lower"]
      ),
    
    mis_sek_k_upper =
      unname(
        sek_selected_k_range["upper"]
      ),
    
    mis_sek_selected_sgn =
      sek_selected_sgn,
    
    mis_sek_support_complete =
      sek_support_complete,
    
    mis_sek_direction_stable =
      sek_direction_stable,
    
    mis_sek_denominator_safe =
      sek_denominator_safe,
    
    mis_sek_formal_guarantee_eligible =
      sek_formal_guarantee_eligible,
    
    mis_sek_instability_reasons =
      sek_instability_reasons,
    
    mis_sek_error =
      mis_sek_error,
    
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
    coef_mis_peel =
      unname(res_peel_v2["coef"]),
    
    coef_mis_sap =
      unname(res_mis_sap["coef"]),
    
    coef_mis_sek =
      unname(res_mis_sek["coef"]),
    
    coef_mm =
      unname(res_mm["coef"]),
    coef_lts = unname(res_lts["coef"]),
    
    se_full = unname(res_full["se"]),
    se_cd = unname(res_cd["se"]),
    se_lev = unname(res_lev["se"]),
    se_dfb = unname(res_dfb["se"]),
    se_mis_alpha = unname(mis_alpha$result["se"]),
    se_mis_oracle = unname(mis_oracle$result["se"]),
    se_mis_peel = unname(res_peel_v2["se"]),
    se_mis_sap = unname(res_mis_sap["se"]),
    se_mis_sek = unname(res_mis_sek["se"]),
    se_mm = unname(res_mm["se"]),
    se_lts = unname(res_lts["se"]),
    
    bias_full = unname(abs(res_full["coef"] - true_b)),
    bias_cd = unname(abs(res_cd["coef"] - true_b)),
    bias_lev = unname(abs(res_lev["coef"] - true_b)),
    bias_dfb = unname(abs(res_dfb["coef"] - true_b)),
    bias_mis_alpha = unname(abs(mis_alpha$result["coef"] - true_b)),
    bias_mis_oracle = unname(abs(mis_oracle$result["coef"] - true_b)),
    bias_mis_peel =
      unname(
        abs(res_peel_v2["coef"] - true_b)
      ),
    
    bias_mis_sap =
      unname(
        abs(res_mis_sap["coef"] - true_b)
      ),
    
    bias_mis_sek =
      unname(
        abs(res_mis_sek["coef"] - true_b)
      ),
    bias_mm = unname(abs(res_mm["coef"] - true_b)),
    bias_lts = unname(
      abs(res_lts["coef"] - true_b)
    ),
    
    deletion_interval_type =
      "naive_post_deletion_ols",
    
    cov_full = check_coverage_v2(
      res_full["coef"],
      res_full["se"],
      true_b
    ),
    
    cov_cd = check_coverage_v2(
      res_cd["coef"],
      res_cd["se"],
      true_b
    ),
    
    cov_lev = check_coverage_v2(
      res_lev["coef"],
      res_lev["se"],
      true_b
    ),
    
    cov_dfb = check_coverage_v2(
      res_dfb["coef"],
      res_dfb["se"],
      true_b
    ),
    
    cov_mis_alpha = check_coverage_v2(
      mis_alpha$result["coef"],
      mis_alpha$result["se"],
      true_b
    ),
    
    cov_mis_oracle = check_coverage_v2(
      mis_oracle$result["coef"],
      mis_oracle$result["se"],
      true_b
    ),
    
    cov_mis_peel = check_coverage_v2(
      res_peel_v2["coef"],
      res_peel_v2["se"],
      true_b
    ),
    
    cov_mis_sap = check_coverage_v2(
      res_mis_sap["coef"],
      res_mis_sap["se"],
      true_b
    ),
    
    cov_mis_sek = check_coverage_v2(
      res_mis_sek["coef"],
      res_mis_sek["se"],
      true_b
    ),
    
    cov_mm = check_coverage_v2(
      res_mm["coef"],
      res_mm["se"],
      true_b
    ),
    
    cov_lts = check_coverage_v2(
      res_lts["coef"],
      res_lts["se"],
      true_b
    ),
    
    cpu_full = unname(cpu_full),
    cpu_cd = unname(cpu_cd),
    cpu_lev = unname(cpu_lev),
    cpu_dfb = unname(cpu_dfb),
    cpu_mis_alpha = unname(cpu_mis_alpha),
    cpu_mis_oracle = unname(cpu_mis_oracle),
    cpu_peel_v2 = unname(cpu_peel_v2),
    cpu_mis_sap =
      unname(cpu_mis_sap),
    
    cpu_mis_sek =
      unname(cpu_mis_sek),
    cpu_mm = unname(cpu_mm),
    cpu_lts = unname(cpu_lts),
    
    stringsAsFactors = FALSE
  )
}
