# ==============================================================================
# File: /R/sim_robust_engine.R
# Purpose: Core iteration engine for the robust comparison simulation. Generates
#          one Monte Carlo draw: creates clean data, injects contamination,
#          fits all 7 estimators (OLS, CD-cleaned, Leverage-cleaned, DFBETAS-cleaned,
#          MIS-cleaned, MM, LTS), and returns a flat 1-row data.frame of bias,
#          raw coefficients, and 95% CI coverage for robust row-binding across
#          parallel simulation loops.
# ==============================================================================

#' Check 95% Wald Interval Coverage
#'
#' @param coef Numeric; the point estimate of the coefficient.
#' @param se Numeric; the standard error of the coefficient estimate.
#' @param true_b Numeric; the true population parameter value.
#'
#' @return Integer: 1L if the true value falls within the 95% CI, 0L if not,
#'         NA_integer_ if either input is NA.
check_coverage <- function(coef, se, true_b) {
  if (is.na(coef) || is.na(se)) return(NA_integer_)
  lo <- coef - 1.96 * se
  hi <- coef + 1.96 * se
  as.integer(true_b >= lo & true_b <= hi)
}

#' Single Iteration of the Robust Comparison Simulation
#'
#' @param iter Integer; the current iteration index.
#' @param n Integer; sample size (default = 1000).
#' @param p Integer; number of predictors (default = 1).
#' @param x_type Character; distribution for the design matrix X. One of
#'        "normal", "mixed_normal", "contaminated".
#' @param error_type Character; distribution for the error term. One of
#'        "normal", "mixed_normal", "skewed_t", "golm", "beta_logistic",
#'        "gpd", "contaminated", "pareto".
#' @param outlier_method Character; contamination topology. One of "none",
#'        "vertical_outlier", "good_leverage", "bad_leverage".
#' @param k Integer; number of observations to contaminate.
#' @param magnitude Numeric; severity multiplier for the injected shift.
#'
#' @return A 1-row data.frame containing iteration metadata, absolute bias,
#'         raw coefficients, and CI coverage flags for all 7 estimators.
run_robust_comparison_iter <- function(iter,
                                       n = 1000,
                                       p = 1,
                                       x_type = "normal",
                                       error_type = "normal",
                                       outlier_method,
                                       k,
                                       magnitude) {
  
  # ---------------------------------------------------------
  # 1. Data Generation & Injection
  # ---------------------------------------------------------
  dat_clean <- generate_complex_data(
    n = n, p = p,
    x_type = x_type,
    error_type = error_type
  )
  true_b <- dat_clean$true_beta[1]
  
  if (outlier_method != "none") {
    dat <- apply_influence_shift(
      dat_clean,
      method = outlier_method,
      k = k,
      magnitude = magnitude
    )
  } else {
    dat <- dat_clean
  }
  df <- data.frame(y = dat$y, x = dat$X[, 1])
  
  # ---------------------------------------------------------
  # 2. Baseline Contaminated OLS
  # ---------------------------------------------------------
  res_full <- fit_clean_ols(y ~ x, data = df, exclude_idx = integer(0))
  mod_full <- stats::lm(y ~ x, data = df)
  beta_full <- unname(stats::coef(mod_full)["x"])
  
  # ---------------------------------------------------------
  # 3. Classical Diagnostics (Using Default Statistical Thresholds)
  # Passing k = NULL forces get_classical_set to use theoretical defaults
  # (4/n for Cook's D, 2p/n for Leverage, 2/sqrt(n) for DFBETAS)
  # ---------------------------------------------------------
  cd_idx  <- get_classical_set(mod_full, target_var = "x", k = NULL, metric = "cooks_d")
  lev_idx <- get_classical_set(mod_full, target_var = "x", k = NULL, metric = "leverage")
  dfb_idx <- get_classical_set(mod_full, target_var = "x", k = NULL, metric = "dfbetas_target")
  
  res_cd  <- fit_clean_ols(y ~ x, data = df, exclude_idx = cd_idx)
  res_lev <- fit_clean_ols(y ~ x, data = df, exclude_idx = lev_idx)
  res_dfb <- fit_clean_ols(y ~ x, data = df, exclude_idx = dfb_idx)
  
  # ---------------------------------------------------------
  # 4. Direct Robust Estimation (MM and LTS)
  # ---------------------------------------------------------
  res_mm  <- fit_mm_estimator(y ~ x, data = df)
  res_lts <- fit_lts_estimator(y ~ x, data = df)
  
  # ---------------------------------------------------------
  # 5. MIS: Determine Dynamic 'k' via Robust Scale
  # We use the MM-estimator strictly to *count* the number of
  # severe structural outliers, rather than looking at its slope.
  # ---------------------------------------------------------
  mod_mm_obj <- tryCatch(
    robustbase::lmrob(y ~ x, data = df, setting = "KS2014"),
    error = function(e) NULL
  )
  
  if (!is.null(mod_mm_obj)) {
    robust_std_res <- mod_mm_obj$residuals / mod_mm_obj$scale
    dynamic_k <- sum(abs(robust_std_res) > 3, na.rm = TRUE)
  } else {
    dynamic_k <- 0
  }
  
  # ---------------------------------------------------------
  # 6. MIS: Direction Selection via Max Coefficient Shift
  # ---------------------------------------------------------
  if (dynamic_k == 0) {
    # CLEAN DATA: No points removed.
    res_mis <- res_full
  } else {
    # Extract the top-k indices for both directions
    mis_idx_pos <- fast_sens_topk(mod_full, pos = 2, sign =  1, k = dynamic_k)
    mis_idx_neg <- fast_sens_topk(mod_full, pos = 2, sign = -1, k = dynamic_k)
    
    # Fit cleaned OLS for both directions
    res_pos <- fit_clean_ols(y ~ x, data = df, exclude_idx = mis_idx_pos)
    res_neg <- fit_clean_ols(y ~ x, data = df, exclude_idx = mis_idx_neg)
    
    # Compute absolute coefficient shift from the contaminated baseline
    shift_pos <- abs(res_pos["coef"] - beta_full)
    shift_neg <- abs(res_neg["coef"] - beta_full)
    
    # Handle NA safely: if one direction failed, take the other
    if (is.na(shift_pos) && is.na(shift_neg)) {
      res_mis <- res_full
    } else if (is.na(shift_pos)) {
      res_mis <- res_neg
    } else if (is.na(shift_neg)) {
      res_mis <- res_pos
    } else {
      # Pick the direction whose removal shifts the coefficient the most
      res_mis <- if (shift_pos >= shift_neg) res_pos else res_neg
    }
  }
  
  # ---------------------------------------------------------
  # 7. Compile Bias and Coverage Metrics
  # ---------------------------------------------------------
  res <- data.frame(
    iter           = iter,
    x_type         = x_type,
    error_type     = error_type,
    outlier_method = outlier_method,
    set_size       = if (outlier_method == "none") 0L else k,
    dynamic_k      = dynamic_k,
    
    # Absolute Bias
    bias_full = unname(abs(res_full["coef"] - true_b)),
    bias_cd   = unname(abs(res_cd["coef"]   - true_b)),
    bias_lev  = unname(abs(res_lev["coef"]  - true_b)),
    bias_dfb  = unname(abs(res_dfb["coef"]  - true_b)),
    bias_mis  = unname(abs(res_mis["coef"]  - true_b)),
    bias_mm   = unname(abs(res_mm["coef"]   - true_b)),
    bias_lts  = unname(abs(res_lts["coef"]  - true_b)),
    
    # 95% CI Coverage
    cov_full = check_coverage(res_full["coef"], res_full["se"], true_b),
    cov_cd   = check_coverage(res_cd["coef"],   res_cd["se"],   true_b),
    cov_lev  = check_coverage(res_lev["coef"],  res_lev["se"],  true_b),
    cov_dfb  = check_coverage(res_dfb["coef"],  res_dfb["se"],  true_b),
    cov_mis  = check_coverage(res_mis["coef"],  res_mis["se"],  true_b),
    cov_mm   = check_coverage(res_mm["coef"],   res_mm["se"],   true_b),
    cov_lts  = check_coverage(res_lts["coef"],  res_lts["se"],  true_b),
    
    # Raw Coefficients
    coef_full = unname(res_full["coef"]),
    coef_cd   = unname(res_cd["coef"]),
    coef_lev  = unname(res_lev["coef"]),
    coef_dfb  = unname(res_dfb["coef"]),
    coef_mis  = unname(res_mis["coef"]),
    coef_mm   = unname(res_mm["coef"]),
    coef_lts  = unname(res_lts["coef"]),
    
    stringsAsFactors = FALSE
  )
  
  return(res)
}