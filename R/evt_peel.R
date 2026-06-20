# ==============================================================================
# File: /R/evt_peel.R
# Purpose: Purpose-built EVT test for the iterative peel. Combines all four
#          components needed for correct sequential testing of detected sets:
#
#          1. Exact Dinkelbach block maxima  (from evt_iter_dm)
#          2. Marginal tail analysis         (from Paper 1, §3.3)
#          3. Robust GEV fitting cascade     (from fit_gev_robust)
#          4. M-scaling                      (from evt_iter, missing in evt_iter_dm)
#
# Why a separate function:
#   evt_iter_dm is missing M-scaling, causing catastrophic over-rejection
#   when testing the detected MIS (selection bias). evt_iter has M-scaling
#   but uses greedy block maxima. evt_iter_dm_v2 "constrained" has the
#   tail analysis but no M-scaling. No single existing function combines
#   all four components. This function does.
#
# Dependencies: helpers_local.R (fwl, make_blocks, dfbeta_numeric),
#               exact_dfb_bmx.R (exact_dfb_bmx),
#               evt_iter_dm.R (fit_gev_robust)
# ==============================================================================


#' EVT Test for Iterative Peel (Full Pipeline)
#'
#' Tests whether a given set's influence on a target coefficient is
#' statistically excessive, using the complete Paper 1 pipeline with
#' exact block maxima and correct M-scaling.
#'
#' Pipeline:
#'   1. FWL orthogonalization → isolate target coefficient
#'   2. Compute true DFBETA of the test set
#'   3. Exact block maxima via Dinkelbach on non-set observations
#'   4. Marginal tail analysis: fit GEV to block max of |X_fwl| and |R_fwl|
#'      separately, test for Fréchet domain, determine tail coefficient
#'   5. Fit GEV to |DFBETA block maxima| with constrained shape = tail_coef;
#'      fall back to robust cascade if constrained fit fails
#'   6. M-scale: transform single-block GEV → max-of-M-blocks GEV
#'   7. Compute p-value against the M-scaled distribution
#'
#' @param y Numeric vector; response variable.
#' @param x Numeric vector; target predictor (the one whose coefficient
#'        we're testing influence on).
#' @param Z Numeric matrix; all other regressors (intercept + controls).
#'        Must NOT contain x.
#' @param set Integer vector; indices of the set being tested.
#' @param block_count Integer; number of blocks for block maxima.
#'        Default = 20.
#'
#' @return A 1-row data.frame with columns:
#'   \item{shape}{Numeric; GEV shape parameter (before M-scaling; shape
#'         is invariant to the M-transformation).}
#'   \item{scale}{Numeric; GEV scale parameter (before M-scaling).}
#'   \item{loc}{Numeric; GEV location parameter (before M-scaling).}
#'   \item{set_dfb}{Numeric; the observed DFBETA of the test set.}
#'   \item{p_value}{Numeric; p-value against the M-scaled GEV. This is
#'         the probability that the maximum of M block maxima exceeds
#'         the observed |set_dfb|.}
#'   \item{converged}{Logical; TRUE if the GEV fit succeeded and p_value
#'         is finite.}
#'   \item{tail_coef}{Numeric; the constrained shape parameter from
#'         marginal tail analysis (0 if Gumbel domain).}
#'   \item{tail_source}{Character; which marginal determined the tail:
#'         "X", "R", "both", or "gumbel".}
#' @export
evt_peel <- function(y, x, Z, set, block_count = 20L) {
  
  # Failure template
  fail_row <- data.frame(
    shape       = NA_real_,
    scale       = NA_real_,
    loc         = NA_real_,
    set_dfb     = NA_real_,
    p_value     = NA_real_,
    converged   = FALSE,
    tail_coef   = NA_real_,
    tail_source = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # ==================================================================
  # 1. FWL Orthogonalization
  # ==================================================================
  fwl_vars <- fwl(y = y, X = x, Z = Z)
  Y_fwl <- fwl_vars[, 1]
  X_fwl <- fwl_vars[, 2]
  R_fwl <- residuals(lm(Y_fwl ~ X_fwl - 1))
  
  # ==================================================================
  # 2. True DFBETA of the test set
  # ==================================================================
  set_dfb <- dfbeta_numeric(Y_fwl, cbind(X_fwl), set, col_X = 1L)
  
  # ==================================================================
  # 3. Exact block maxima via Dinkelbach
  # ==================================================================
  bmx <- tryCatch(
    exact_dfb_bmx(X = X_fwl, R = R_fwl, set = set,
                  block_count = block_count),
    error = function(e) NULL
  )
  if (is.null(bmx) || length(bmx) < 3L) return(fail_row)
  
  M <- length(bmx)  # actual number of blocks (may differ from requested)
  
  # ==================================================================
  # 4. Marginal tail analysis (Paper 1, Section 3.3)
  #
  #    Estimate the tail behavior of X_fwl and R_fwl separately.
  #    If either is in the Fréchet domain (ξ > 0 with statistical
  #    significance), the DFBETA = X·R products inherit that tail,
  #    and we constrain the DFBETA GEV shape accordingly.
  # ==================================================================
  X_nonsig <- X_fwl[-set]
  R_nonsig <- R_fwl[-set]
  block_size_marginal <- length(X_nonsig) %/% block_count
  
  # Default: Gumbel domain (light tails)
  tail_coef   <- 0
  tail_source <- "gumbel"
  
  if (block_size_marginal >= 2L) {
    # Block maxima of |X| and |R| for marginal tail estimation
    bm_X <- tryCatch({
      apply(make_blocks(abs(X_nonsig), block_size_marginal), 2, max)
    }, error = function(e) NULL)
    
    bm_R <- tryCatch({
      apply(make_blocks(abs(R_nonsig), block_size_marginal), 2, max)
    }, error = function(e) NULL)
    
    # Fit marginal GEVs
    xi_x <- 0
    xi_x_sig <- FALSE
    if (!is.null(bm_X) && length(bm_X) >= 3L) {
      x_evd <- tryCatch(evd::fgev(bm_X), error = function(e) NULL)
      if (!is.null(x_evd)) {
        xi_x <- x_evd$estimate["shape"]
        # Fréchet test: lower 95% CI bound > 0
        xi_x_sig <- (xi_x - 1.96 * x_evd$std.err["shape"]) > 0
      }
    }
    
    xi_r <- 0
    xi_r_sig <- FALSE
    if (!is.null(bm_R) && length(bm_R) >= 3L) {
      r_evd <- tryCatch(evd::fgev(bm_R), error = function(e) NULL)
      if (!is.null(r_evd)) {
        xi_r <- r_evd$estimate["shape"]
        xi_r_sig <- (xi_r - 1.96 * r_evd$std.err["shape"]) > 0
      }
    }
    
    # Determine tail coefficient
    if (xi_x_sig && xi_r_sig) {
      tail_coef   <- max(xi_x, xi_r)
      tail_source <- "both"
    } else if (xi_x_sig) {
      tail_coef   <- xi_x
      tail_source <- "X"
    } else if (xi_r_sig) {
      tail_coef   <- xi_r
      tail_source <- "R"
    }
    # else: both Gumbel → tail_coef stays 0, tail_source stays "gumbel"
  }
  
  # ==================================================================
  # 5. Fit GEV to |DFBETA block maxima|
  #
  #    Primary: constrained shape = tail_coef (Paper 1 approach)
  #    Fallback: robust cascade (if constrained fit fails)
  # ==================================================================
  abs_bmx <- abs(bmx)
  fit_evd <- NULL
  
  # 5a. Try constrained fit (shape fixed at tail_coef)
  fit_evd <- tryCatch(
    evd::fgev(abs_bmx, shape = tail_coef),
    error = function(e) NULL
  )
  
  # Validate: scale must be positive
  if (!is.null(fit_evd) && fit_evd$estimate["scale"] <= 0) {
    fit_evd <- NULL
  }
  
  # 5b. Fallback: robust cascade (unconstrained, from fit_gev_robust)
  if (is.null(fit_evd)) {
    fit_evd <- tryCatch(
      fit_gev_robust(abs_bmx),
      error = function(e) NULL
    )
  }
  
  # 5c. Both failed
  if (is.null(fit_evd) || fit_evd$estimate["scale"] <= 0) {
    return(fail_row)
  }
  
  sigma <- unname(fit_evd$estimate["scale"])
  mu    <- unname(fit_evd$estimate["loc"])
  xi    <- if ("shape" %in% names(fit_evd$estimate)) {
    unname(fit_evd$estimate["shape"])
  } else {
    tail_coef
  }
  
  # ==================================================================
  # 6. M-Scaling: single-block GEV → max-of-M-blocks GEV
  #
  #    If X₁,...,X_M ~ iid GEV(μ, σ, ξ), then
  #      max(X₁,...,X_M) ~ GEV(μ_M, σ_M, ξ)
  #
  #    Fréchet (ξ > 0):  μ_M = μ + σ(M^ξ − 1)/ξ,   σ_M = σ·M^ξ
  #    Gumbel  (ξ ≈ 0):  μ_M = μ + σ·log(M),        σ_M = σ
  #    Weibull (ξ < 0):  μ_M = μ + σ(M^ξ − 1)/ξ,   σ_M = σ·M^ξ
  #                       (same formula, M^ξ < 1 so σ shrinks)
  # ==================================================================
  if (abs(xi) > 1e-6) {
    mu_M    <- mu + sigma * (M^xi - 1) / xi
    sigma_M <- sigma * M^xi
  } else {
    mu_M    <- mu + sigma * log(M)
    sigma_M <- sigma
  }
  
  # Guard: sigma_M must remain positive
  if (sigma_M <= 0) return(fail_row)
  
  # ==================================================================
  # 7. P-value against the M-scaled distribution
  # ==================================================================
  p_val <- 1 - evd::pgev(
    q     = abs(set_dfb),
    loc   = mu_M,
    scale = sigma_M,
    shape = xi
  )
  if (!is.finite(p_val)) p_val <- NA_real_
  
  # ==================================================================
  # 8. Return
  # ==================================================================
  data.frame(
    shape       = xi,
    scale       = sigma,
    loc         = mu,
    set_dfb     = unname(set_dfb),
    p_value     = unname(p_val),
    converged   = is.finite(p_val),
    tail_coef   = tail_coef,
    tail_source = tail_source,
    stringsAsFactors = FALSE
  )
}