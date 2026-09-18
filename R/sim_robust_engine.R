# ==============================================================================
# File: /R/sim_robust_engine.R
# Purpose: Core iteration engine for the robust comparison simulation. Generates
#          one Monte Carlo draw: creates clean data, injects contamination,
#          fits all estimators, and returns a flat 1-row data.frame of bias,
#          raw coefficients, and 95% CI coverage for robust row-binding across
#          parallel simulation loops.
#
# Final Script 04 comparison methods:
#   1. OLS
#   2. Leverage deletion + OLS refit
#   3. Cook's-distance deletion + OLS refit
#   4. DFBETAS deletion + OLS refit
#   5. MIS with oracle k and oracle direction + OLS refit
#   6. MM
#   7. LTS
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
  if (
    !is.finite(coef) ||
    !is.finite(se) ||
    se < 0
  ) {
    return(NA_integer_)
  }
  lo <- coef - 1.96 * se
  hi <- coef + 1.96 * se
  as.integer(true_b >= lo & true_b <= hi)
}

compute_overlap <- function(detected, true_idx) {
  if (is.null(true_idx) || length(true_idx) == 0L) {
    return(NA_real_)
  }
  
  detected <- detected[is.finite(detected)]
  
  if (length(detected) == 0L) {
    return(0)
  }
  
  length(intersect(detected, true_idx)) / length(true_idx)
}

extract_target_lm <- function(
    model,
    target_var = "x"
) {
  tryCatch(
    {
      sm <- summary(model)$coefficients
      
      if (!target_var %in% rownames(sm)) {
        return(
          c(
            coef = NA_real_,
            se = NA_real_
          )
        )
      }
      
      c(
        coef = unname(
          sm[target_var, "Estimate"]
        ),
        se = unname(
          sm[target_var, "Std. Error"]
        )
      )
    },
    error = function(e) {
      c(
        coef = NA_real_,
        se = NA_real_
      )
    }
  )
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
#'         raw coefficients, and CI coverage flags for all estimators.
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
    true_idx <- dat$outlier_indices
  } else {
    dat <- dat_clean
    true_idx <- integer(0)
  }
  df <- data.frame(y = dat$y, x = dat$X[, 1])
  
  # ---------------------------------------------------------
  # 2. OLS
  # ---------------------------------------------------------
  t0 <- proc.time()[3L]
  
  mod_full <- stats::lm(
    y ~ x,
    data = df
  )
  
  res_full <- extract_target_lm(
    mod_full,
    target_var = "x"
  )
  
  cpu_full <- proc.time()[3L] - t0
  
  # ---------------------------------------------------------
  # 3. Shared model quantities
  # ---------------------------------------------------------
  
  n_model <- stats::nobs(mod_full)
  p_model <- length(stats::coef(mod_full))
  
  
  # ---------------------------------------------------------
  # 4. Leverage deletion + OLS refit
  # ---------------------------------------------------------
  
  t0 <- proc.time()[3L]
  
  lev_values <- get_leverage(
    mod_full
  )
  
  lev_idx <- which(
    is.finite(lev_values) &
      lev_values > (2 * p_model / n_model)
  )
  
  res_lev <- fit_clean_ols(
    y ~ x,
    data = df,
    exclude_idx = lev_idx
  )
  
  cpu_lev <- unname(
    proc.time()[3L] - t0
  )
  
  
  # ---------------------------------------------------------
  # 5. Cook's-distance deletion + OLS refit
  # ---------------------------------------------------------
  
  t0 <- proc.time()[3L]
  
  cd_values <- get_cooks_d(
    mod_full
  )
  
  cd_idx <- which(
    is.finite(cd_values) &
      cd_values > (4 / n_model)
  )
  
  res_cd <- fit_clean_ols(
    y ~ x,
    data = df,
    exclude_idx = cd_idx
  )
  
  cpu_cd <- unname(
    proc.time()[3L] - t0
  )
  
  
  # ---------------------------------------------------------
  # 6. DFBETAS deletion + OLS refit
  # ---------------------------------------------------------
  
  t0 <- proc.time()[3L]
  
  dfb_values <- get_dfbetas(
    mod_full,
    target_var = "x"
  )
  
  dfb_idx <- which(
    is.finite(dfb_values) &
      abs(dfb_values) > (2 / sqrt(n_model))
  )
  
  res_dfb <- fit_clean_ols(
    y ~ x,
    data = df,
    exclude_idx = dfb_idx
  )
  
  cpu_dfb <- unname(
    proc.time()[3L] - t0
  )
  
  # ---------------------------------------------------------
  # 7. MIS with oracle k and oracle direction
  # ---------------------------------------------------------
  t0 <- proc.time()[3L]
  
  if (outlier_method == "none") {
    
    k_oracle_val <- 0L
    mis_oracle_idx <- integer(0)
    oracle_dfbeta_delta <- 0
    oracle_sign <- 0L
    res_mis_oracle <- res_full
    
  } else {
    
    k_oracle_val <- as.integer(k)
    
    oracle_clean <- fit_clean_ols(
      y ~ x,
      data = df,
      exclude_idx = true_idx
    )
    
    # Dinkelbach's directional objective uses the deletion-DFBETA
    # orientation:
    #
    #   beta_full - beta_deleted
    #
    # Therefore the oracle direction must use the same orientation.
    oracle_dfbeta_delta <- unname(
      res_full["coef"] - oracle_clean["coef"]
    )
    
    if (!is.finite(oracle_dfbeta_delta)) {
      stop(
        "Oracle direction could not be determined: ",
        "non-finite oracle DFBETA difference."
      )
    }
    
    oracle_sign <- if (
      oracle_dfbeta_delta >= 0
    ) {
      1L
    } else {
      -1L
    }
    
    mis_oracle_idx <- dinkelbach_topk_lm(
      mod = mod_full,
      pos = 2L,
      sign = oracle_sign,
      k = k_oracle_val
    )
    
    res_mis_oracle <- fit_clean_ols(
      y ~ x,
      data = df,
      exclude_idx = mis_oracle_idx
    )
  }
  
  cpu_mis_oracle <- proc.time()[3L] - t0
  
  
  # ---------------------------------------------------------
  # 8. MM
  # ---------------------------------------------------------
  t0 <- proc.time()[3L]
  
  res_mm <- fit_mm_estimator(
    y ~ x,
    data = df
  )
  
  cpu_mm <- proc.time()[3L] - t0
  
  
  # ---------------------------------------------------------
  # 9. LTS
  # ---------------------------------------------------------
  t0 <- proc.time()[3L]
  
  res_lts <- fit_lts_estimator(
    y ~ x,
    data = df
  )
  
  cpu_lts <- proc.time()[3L] - t0
  
  
  # ---------------------------------------------------------
  # 10. Injected-set recovery
  # ---------------------------------------------------------
  overlap_lev <- compute_overlap(lev_idx, true_idx)
  overlap_cd <- compute_overlap(cd_idx, true_idx)
  overlap_dfb <- compute_overlap(dfb_idx, true_idx)
  overlap_mis_oracle <- compute_overlap(
    mis_oracle_idx,
    true_idx
  )
  
  # ---------------------------------------------------------
  # 11. Compile Bias and Coverage Metrics
  # ---------------------------------------------------------
  res <- data.frame(
    iter = iter,
    
    n_obs = as.integer(n),
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
    
    # Selected-set sizes
    k_lev = length(lev_idx),
    k_cd = length(cd_idx),
    k_dfb = length(dfb_idx),
    k_oracle = k_oracle_val,
    
    # Injected-set recovery
    overlap_lev = overlap_lev,
    overlap_cd = overlap_cd,
    overlap_dfb = overlap_dfb,
    overlap_mis_oracle = overlap_mis_oracle,
    oracle_direction = oracle_sign,
    oracle_dfbeta_delta = oracle_dfbeta_delta,
    
    deletion_interval_type =
      "naive_post_deletion_ols",
    
    # Coefficients
    coef_full = unname(res_full["coef"]),
    coef_lev = unname(res_lev["coef"]),
    coef_cd = unname(res_cd["coef"]),
    coef_dfb = unname(res_dfb["coef"]),
    coef_mis_oracle =
      unname(res_mis_oracle["coef"]),
    coef_mm = unname(res_mm["coef"]),
    coef_lts = unname(res_lts["coef"]),
    
    # Standard errors
    se_full = unname(res_full["se"]),
    se_lev = unname(res_lev["se"]),
    se_cd = unname(res_cd["se"]),
    se_dfb = unname(res_dfb["se"]),
    se_mis_oracle =
      unname(res_mis_oracle["se"]),
    se_mm = unname(res_mm["se"]),
    se_lts = unname(res_lts["se"]),
    
    # Absolute bias
    bias_full =
      abs(unname(res_full["coef"]) - true_b),
    bias_lev =
      abs(unname(res_lev["coef"]) - true_b),
    bias_cd =
      abs(unname(res_cd["coef"]) - true_b),
    bias_dfb =
      abs(unname(res_dfb["coef"]) - true_b),
    bias_mis_oracle =
      abs(unname(res_mis_oracle["coef"]) - true_b),
    bias_mm =
      abs(unname(res_mm["coef"]) - true_b),
    bias_lts =
      abs(unname(res_lts["coef"]) - true_b),
    
    # 95% coverage
    cov_full = check_coverage(
      res_full["coef"], res_full["se"], true_b
    ),
    cov_lev = check_coverage(
      res_lev["coef"], res_lev["se"], true_b
    ),
    cov_cd = check_coverage(
      res_cd["coef"], res_cd["se"], true_b
    ),
    cov_dfb = check_coverage(
      res_dfb["coef"], res_dfb["se"], true_b
    ),
    cov_mis_oracle = check_coverage(
      res_mis_oracle["coef"],
      res_mis_oracle["se"],
      true_b
    ),
    cov_mm = check_coverage(
      res_mm["coef"], res_mm["se"], true_b
    ),
    cov_lts = check_coverage(
      res_lts["coef"], res_lts["se"], true_b
    ),
    
    # Runtime
    cpu_full = cpu_full,
    cpu_lev = cpu_lev,
    cpu_cd = cpu_cd,
    cpu_dfb = cpu_dfb,
    cpu_mis_oracle = cpu_mis_oracle,
    cpu_mm = cpu_mm,
    cpu_lts = cpu_lts,
    
    stringsAsFactors = FALSE
  )
  
  return(res)
}