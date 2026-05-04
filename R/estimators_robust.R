# ==============================================================================
# File: /R/estimators_robust.R
# Purpose: Wrappers for robust regression estimators (LTS and MM) and a clean 
#          OLS fallback for direct cross-comparison in simulations.
# ==============================================================================

# Ensure the required package is available without strictly attaching it to the search path
if (!requireNamespace("robustbase", quietly = TRUE)) {
  stop("The 'robustbase' package is required. Please install it using install.packages('robustbase').")
}

#' Fit Robust MM-Estimator safely
#'
#' @param formula An object of class "formula"
#' @param data A data frame containing the variables in the model
#' @return A numeric vector of coefficients, or properly sized NAs if convergence fails
fit_mm_estimator <- function(formula, data) {
  tryCatch({
    mod <- robustbase::lmrob(formula, data = data, setting = "KS2014")
    smry <- summary(mod)$coefficients
    # Return coef and se for the slope 'x', matching fit_clean_ols output format
    if ("x" %in% rownames(smry)) {
      c(coef = smry["x", "Estimate"], se = smry["x", "Std. Error"])
    } else {
      c(coef = NA_real_, se = NA_real_)
    }
  }, error = function(e) {
    warning("MM-estimator failed to converge: ", conditionMessage(e))
    c(coef = NA_real_, se = NA_real_)
  })
}

#' Fit Least Trimmed Squares (LTS) Estimator safely
#'
#' @param formula An object of class "formula"
#' @param data A data frame containing the variables in the model
#' @return A numeric vector of coefficients, or properly sized NAs if convergence fails
fit_lts_estimator <- function(formula, data) {
  tryCatch({
    # ltsReg computes the Fast-LTS estimator. 
    mod <- robustbase::ltsReg(formula, data = data)
    coefs <- stats::coef(mod)
    # LTS does not provide classical SEs; compute them from the reweighted
    # clean subset OLS (observations with weight == 1 in the final step)
    se_val <- tryCatch({
      # ltsReg stores raw residuals; the reweighted fit uses non-outlying obs
      raw_res <- mod$residuals
      wts <- mod$lts.wt  # binary weights: 1 = clean, 0 = outlier
      if (!is.null(wts) && sum(wts) > length(coefs)) {
        X_mat <- stats::model.matrix(formula, data)
        X_clean <- X_mat[wts == 1, , drop = FALSE]
        rss_clean <- sum(raw_res[wts == 1]^2)
        n_clean <- sum(wts)
        p_mod <- ncol(X_clean)
        sigma2 <- rss_clean / (n_clean - p_mod)
        se_vec <- sqrt(diag(sigma2 * solve(crossprod(X_clean))))
        names(se_vec) <- colnames(X_clean)
        if ("x" %in% names(se_vec)) se_vec["x"] else NA_real_
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)
    
    # Return coef and se for the slope 'x', matching fit_clean_ols output format
    if ("x" %in% names(coefs)) {
      c(coef = unname(coefs["x"]), se = unname(se_val))
    } else {
      c(coef = NA_real_, se = NA_real_)
    }
  }, error = function(e) {
    warning("LTS estimator failed: ", conditionMessage(e))
    c(coef = NA_real_, se = NA_real_)
  })
}

#' Fit OLS on a dynamically cleaned subset (Returns Coef and SE)
fit_clean_ols <- function(formula, data, exclude_idx) {
  tryCatch({
    # Safely drop NA indices if any leaked in from the combinatorial search
    exclude_idx <- exclude_idx[!is.na(exclude_idx)]
    
    if (length(exclude_idx) > 0) {
      clean_data <- data[-exclude_idx, , drop = FALSE]
    } else {
      clean_data <- data
    }
    
    mod <- stats::lm(formula, data = clean_data)
    smry <- summary(mod)$coefficients
    
    # Check if 'x' survived the fit (didn't get dropped for singularity)
    if ("x" %in% rownames(smry)) {
      c(coef = smry["x", "Estimate"], se = smry["x", "Std. Error"])
    } else {
      c(coef = NA_real_, se = NA_real_)
    }
  }, error = function(e) {
    # If the matrix inversion fails entirely, return NAs safely
    c(coef = NA_real_, se = NA_real_)
  })
}