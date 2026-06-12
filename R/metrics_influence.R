# ==============================================================================
# File: /R/metrics_influence.R
# Purpose: Evaluation metrics for MIS detection quality beyond set overlap.
#          Provides influence ratio (magnitude-based), nestedness tracing
#          (set stability across k), null calibration statistics (p-value
#          uniformity), and GEV shape parameter accuracy assessment.
#
# Dependencies: Requires helpers_local.R (provides dfbeta_numeric),
#               dinkelbach_topk.R (provides dinkelbach_topk_lm,
#               dinkelbach_topk_refined), fast_sens_topk.R (provides
#               fast_sens_topk).
# ==============================================================================


#' Influence Ratio Between Detected and True Sets
#'
#' Computes the ratio of absolute DFBETA magnitudes: |DFBETA(detected)| /
#' |DFBETA(true)|. Values near 1.0 mean the detected set captures nearly
#' all the true influence, even if set membership differs. Values above 1.0
#' mean the detected set is MORE influential than the true set (possible
#' when the detection algorithm finds a better set than the injected one).
#'
#' @param y     Numeric vector; FWL-projected response (Y_fwl).
#' @param x_fwl Numeric vector; FWL-projected predictor (X_fwl).
#' @param detected_set Integer vector; indices of the detected set.
#' @param true_set     Integer vector; indices of the true (injected) set.
#'
#' @return Numeric scalar in [0, Inf). Returns NA_real_ if the true set
#'         DFBETA is zero (no influence to compare against).
#' @export
influence_ratio_fwl <- function(y, x_fwl, detected_set, true_set) {
  
  dfb_detected <- dfbeta_numeric(y, cbind(x_fwl), detected_set, col_X = 1L)
  dfb_true     <- dfbeta_numeric(y, cbind(x_fwl), true_set,     col_X = 1L)
  
  if (abs(dfb_true) < 1e-15) return(NA_real_)
  
  abs(dfb_detected) / abs(dfb_true)
}


#' Nestedness Trace Across k = 1, ..., k_max
#'
#' For each k, detects the k-MIS using the specified method and records
#' set composition, Jaccard similarity with the (k-1)-MIS, nestedness
#' violations (S_{k-1} not a subset of S_k), and influence magnitude.
#'
#' @param mod    A fitted lm object.
#' @param pos    Integer; target coefficient position (default = 2).
#' @param k_max  Integer; maximum set size to trace.
#' @param method Character; detection method. One of \code{"dinkelbach"}
#'        (calls \code{dinkelbach_topk_lm}), \code{"greedy"} (calls
#'        \code{fast_sens_topk}), or \code{"dinkelbach_refined"} (calls
#'        \code{dinkelbach_topk_refined}).
#' @param sign   Integer; +1 or -1 direction.
#'
#' @return A data.frame with columns: k, jaccard_with_prev, nested,
#'         influence_magnitude, cpu_seconds.
#' @export
nestedness_trace <- function(mod, pos = 2L, k_max, method = "dinkelbach",
                             sign = 1L) {
  
  detect_fn <- switch(method,
                      "dinkelbach"         = function(k) dinkelbach_topk_lm(mod, pos = pos,
                                                                            sign = sign, k = k),
                      "greedy"             = function(k) fast_sens_topk(mod, pos = pos,
                                                                        sign = sign, k = k),
                      "dinkelbach_refined" = function(k) dinkelbach_topk_refined(mod, pos = pos,
                                                                                 sign = sign,
                                                                                 k = k),
                      stop(sprintf("Unknown method: '%s'.", method))
  )
  
  # FWL components for influence magnitude
  X <- stats::model.matrix(mod)
  y <- stats::model.response(stats::model.frame(mod))
  p <- ncol(X)
  Z_cols <- setdiff(seq_len(p), pos)
  
  if (length(Z_cols) == 0L) {
    x_fwl <- X[, pos]
    y_fwl <- y
  } else {
    Z <- X[, Z_cols, drop = FALSE]
    qr_Z <- qr(Z)
    x_fwl <- qr.resid(qr_Z, X[, pos])
    y_fwl <- qr.resid(qr_Z, y)
  }
  
  # Trace
  prev_set <- integer(0)
  results  <- vector("list", k_max)
  
  for (k in seq_len(k_max)) {
    t0      <- proc.time()[3]
    cur_set <- detect_fn(k)
    cpu     <- proc.time()[3] - t0
    
    # Jaccard similarity with previous k
    if (k == 1L) {
      jaccard <- NA_real_
      is_nested <- NA
    } else {
      intersection_size <- length(intersect(prev_set, cur_set))
      union_size        <- length(union(prev_set, cur_set))
      jaccard    <- intersection_size / union_size
      is_nested  <- all(prev_set %in% cur_set)
    }
    
    # Influence magnitude of the current set
    infl_mag <- abs(dfbeta_numeric(y_fwl, cbind(x_fwl), cur_set, col_X = 1L))
    
    results[[k]] <- data.frame(
      k                   = k,
      jaccard_with_prev   = jaccard,
      nested              = is_nested,
      influence_magnitude = infl_mag,
      cpu_seconds         = cpu,
      stringsAsFactors    = FALSE
    )
    
    prev_set <- cur_set
  }
  
  do.call(rbind, results)
}


#' Null Calibration Statistics for P-Values
#'
#' Evaluates whether a vector of p-values (from null simulations where no
#' outliers are injected) is consistent with Uniform(0,1). This is the
#' fundamental validity check for the EVT test.
#'
#' @param p_values Numeric vector of p-values. NA values are removed.
#'
#' @return A list with components:
#'   \item{ks_pvalue}{Numeric; p-value from a Kolmogorov-Smirnov test
#'                    against Uniform(0,1). Non-significant (> 0.05)
#'                    indicates good calibration.}
#'   \item{rejection_rates}{Named numeric vector; empirical rejection
#'                          rates at alpha = 0.01, 0.05, 0.10.}
#'   \item{qq_data}{Data.frame with columns \code{theoretical} and
#'                  \code{empirical} for QQ plotting.}
#'   \item{n_valid}{Integer; count of non-NA p-values used.}
#' @export
pvalue_calibration_stats <- function(p_values) {
  
  p_clean <- p_values[is.finite(p_values)]
  n_valid <- length(p_clean)
  
  if (n_valid < 3L) {
    return(list(
      ks_pvalue       = NA_real_,
      rejection_rates = c("0.01" = NA_real_, "0.05" = NA_real_,
                          "0.10" = NA_real_),
      qq_data         = data.frame(theoretical = numeric(0),
                                   empirical   = numeric(0)),
      n_valid         = n_valid
    ))
  }
  
  # KS test against Uniform(0,1)
  ks_result <- stats::ks.test(p_clean, "punif", 0, 1)
  
  # Empirical rejection rates
  rej <- c(
    "0.01" = mean(p_clean < 0.01),
    "0.05" = mean(p_clean < 0.05),
    "0.10" = mean(p_clean < 0.10)
  )
  
  # QQ data: sorted empirical vs theoretical quantiles
  p_sorted    <- sort(p_clean)
  theoretical <- (seq_len(n_valid) - 0.5) / n_valid
  
  qq_df <- data.frame(
    theoretical = theoretical,
    empirical   = p_sorted,
    stringsAsFactors = FALSE
  )
  
  list(
    ks_pvalue       = ks_result$p.value,
    rejection_rates = rej,
    qq_data         = qq_df,
    n_valid         = n_valid
  )
}


#' GEV Shape Parameter Accuracy
#'
#' Compares a vector of fitted GEV shape estimates (from Monte Carlo
#' repetitions) against a theoretical prediction. Reports bias, RMSE,
#' and coverage of the theoretical value.
#'
#' @param fitted_shapes   Numeric vector of fitted xi values. NA values
#'                        are removed.
#' @param theoretical_xi  Numeric scalar; the theoretically predicted
#'                        shape parameter (e.g., 0 for Gumbel, or the
#'                        tail index for Frechet).
#'
#' @return A list with components:
#'   \item{mean_bias}{Numeric; mean(fitted - theoretical).}
#'   \item{median_bias}{Numeric; median(fitted - theoretical).}
#'   \item{rmse}{Numeric; root mean squared error.}
#'   \item{coverage_90}{Numeric; proportion of MC draws where the
#'         theoretical value falls within the central 90\% interval
#'         of the fitted distribution (5th to 95th percentile).}
#'   \item{n_valid}{Integer; count of non-NA fitted values used.}
#' @export
shape_parameter_bias <- function(fitted_shapes, theoretical_xi) {
  
  xi_clean <- fitted_shapes[is.finite(fitted_shapes)]
  n_valid  <- length(xi_clean)
  
  if (n_valid < 3L) {
    return(list(
      mean_bias   = NA_real_,
      median_bias = NA_real_,
      rmse        = NA_real_,
      coverage_90 = NA_real_,
      n_valid     = n_valid
    ))
  }
  
  deviations <- xi_clean - theoretical_xi
  
  # Central 90% interval of the fitted distribution
  q05 <- quantile(xi_clean, 0.05)
  q95 <- quantile(xi_clean, 0.95)
  covered <- (theoretical_xi >= q05) && (theoretical_xi <= q95)
  
  list(
    mean_bias   = mean(deviations),
    median_bias = median(deviations),
    rmse        = sqrt(mean(deviations^2)),
    coverage_90 = as.numeric(covered),
    n_valid     = n_valid
  )
}