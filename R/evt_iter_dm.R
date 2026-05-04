# ==============================================================================
# File: /R/evt_iter_dm.R
# Purpose: Wrap the EVT estimation using the exact Dinkelbach's Method (DM) for 
#          block maxima. This bypasses the greedy heuristic in testingMIS to 
#          provide mathematically exact p-values for the Most Influential Set.
#          Implements a multi-attempt GEV fitting strategy (5 attempts):
#            1. Default fgev (BFGS, internal init)
#            2. L-moments starting values (Hosking 1990)
#            3. Conservative Gumbel-like start (ξ ≈ 0)
#            4. Shape-clamped Nelder-Mead (derivative-free, |ξ| ≤ 0.5)
#            5. Median/IQR robust start + Nelder-Mead (50% breakdown point)
# Dependencies: Requires /R/exact_dfb_bmx.R to be sourced.
# ==============================================================================

#' Exact Wrapper for EVT Estimation in Simulations (Dinkelbach's Method)
#'
#' @param y Numeric vector; response variable.
#' @param x Numeric vector; primary predictor variable.
#' @param Z Numeric matrix; covariates to be marginalized out (intercept-only
#'        for simple regression; must NOT contain x).
#' @param set Integer vector; indices of the influential set being tested.
#' @param block_count Integer; number of blocks for block maxima approach 
#'        (default = 20).
#' 
#' @return A 1-row data.frame containing GEV parameters (shape, scale, loc),
#'         the observed set DFBETA, the extreme value p-value, and a 
#'         convergence flag.
#' @importFrom testingMIS dfbeta_numeric
#' @importFrom evd fgev pgev
#' @export
evt_iter_dm <- function(y, x, Z, set, block_count = 20) {
  
  # Failure template — returned when all fitting attempts fail
  fail_row <- data.frame(
    shape   = NA_real_,
    scale   = NA_real_,
    loc     = NA_real_,
    set_dfb = NA_real_,
    p_value = NA_real_,
    converged = FALSE,
    stringsAsFactors = FALSE
  )
  
  # 1. FWL Orthogonalization (Isolating the effect of x on y)
  fwl_vars <- testingMIS:::fwl(y = y, X = x, Z = Z)
  Y_fwl <- fwl_vars[, 1]
  X_fwl <- fwl_vars[, 2]
  
  # Compute the residuals of the orthogonalized model
  R_fwl <- residuals(lm(Y_fwl ~ X_fwl - 1))
  
  # 2. Compute the True DFBETA of the target set
  #    (Pass X_fwl as a 1-column matrix to match dfbeta_numeric's signature)
  set_dfb <- testingMIS::dfbeta_numeric(Y_fwl, cbind(X_fwl), set, col_X = 1L)
  
  # 3. Compute block maxima using exact Dinkelbach's method on actual data
  #    [BUG1 FIX] This replaces the old MC null-draw approach that:
  #    (a) resampled from contaminated residuals (polluted null), and
  #    (b) treated iid rmaxdfbeta draws as block maxima (wrong EVD structure).
  #    Now we use exact_dfb_bmx which divides the non-set observations into
  #    blocks and finds the exact maximum-influence subset within each block
  #    via linear-fractional programming — matching the testingMIS methodology.
  bmx <- tryCatch(
    exact_dfb_bmx(X = X_fwl, R = R_fwl, set = set, block_count = block_count),
    error = function(e) NULL
  )
  if (is.null(bmx) || length(bmx) < 3) return(fail_row)
  
  # 4. Fit GEV to the block maxima
  fit_evd <- tryCatch(fit_gev_robust(abs(bmx)), error = function(e) NULL)
  if (is.null(fit_evd) || fit_evd$estimate["scale"] <= 0) { return(fail_row) }
  
  xi    <- fit_evd$estimate["shape"]
  sigma <- fit_evd$estimate["scale"]
  mu    <- fit_evd$estimate["loc"]
  
  # 5. Compute p-value
  p_val <- 1 - evd::pgev(q = abs(set_dfb), loc = mu, scale = sigma, shape = xi)
  if (!is.finite(p_val)) p_val <- NA_real_
  
  # 6. Return strict 1-row data frame
  data.frame(
    shape     = unname(xi), 
    scale     = unname(sigma),
    loc       = unname(mu),
    set_dfb   = unname(set_dfb), 
    p_value   = unname(p_val), 
    converged = TRUE,
    stringsAsFactors = FALSE
  )
}


#' Robust GEV Fitting with Multiple Fallback Strategies
#'
#' Attempts evd::fgev up to five times with progressively more robust
#' starting values and optimizer switches. Recovers fits that fail under
#' default initialisation due to heavy tails, near-degenerate data, or
#' optimizer sensitivity.
#'
#' Strategy:
#'   1. Default fgev (BFGS, internal starting values)
#'   2. L-moments starting values (Hosking 1990) — robust to heavy tails
#'   3. Conservative Gumbel-like start (ξ ≈ 0)
#'   4. Shape-clamped Nelder-Mead — derivative-free optimizer with ξ
#'      hard-clamped to [-0.5, 0.5] via penalised negative log-likelihood.
#'      Recovers fits where BFGS diverges due to flat or ridged likelihood.
#'   5. Median/IQR robust start with Nelder-Mead — breakdown-resistant
#'      initialisation for heavily contaminated or near-constant block maxima.
#'
#' @param bmx Numeric vector of block maxima (must have length >= 3).
#' @return An fgev fit object, or NULL if all attempts fail.
#' @keywords internal
fit_gev_robust <- function(bmx) {
  
  # Shared validator: scale must be strictly positive
  .valid <- function(fit) {
    !is.null(fit) && fit$estimate["scale"] > 0
  }
  
  # ------------------------------------------------------------------
  # Attempt 1: Default fgev (uses its own MLE starting values)
  # ------------------------------------------------------------------
  fit <- tryCatch(evd::fgev(bmx), error = function(e) NULL)
  if (.valid(fit)) return(fit)
  
  # ------------------------------------------------------------------
  # Attempt 2: L-moments starting values (robust to heavy tails)
  # ------------------------------------------------------------------
  lmom_start <- tryCatch({
    # Simple L-moment estimates for GEV (Hosking 1990)
    n <- length(bmx)
    bmx_sorted <- sort(bmx)
    # L-moment ratios via PWM
    b0 <- mean(bmx_sorted)
    b1 <- sum((seq_len(n) - 1) / (n - 1) * bmx_sorted) / n
    b2 <- sum((seq_len(n) - 1) * (seq_len(n) - 2) / ((n - 1) * (n - 2)) * bmx_sorted) / n
    
    l1 <- b0
    l2 <- 2 * b1 - b0
    t3 <- (6 * b2 - 6 * b1 + b0) / (2 * b1 - b0)
    
    # Approximate ξ from L-skewness (Hosking & Wallis 1997, eq 3.6)
    c_val <- 2 / (3 + t3) - log(2) / log(3)
    xi_est <- 7.8590 * c_val + 2.9554 * c_val^2
    
    if (abs(xi_est) > 0.5) xi_est <- sign(xi_est) * 0.5
    
    gam <- gamma(1 - xi_est)
    sigma_est <- l2 * xi_est / (gam * (2^xi_est - 1))
    mu_est <- l1 - sigma_est * (gam - 1) / xi_est
    
    if (!is.finite(sigma_est) || sigma_est <= 0) {
      sigma_est <- l2 * sqrt(pi) / sqrt(6)
      mu_est <- l1 - 0.5772 * sigma_est
      xi_est <- 0.01
    }
    
    list(loc = mu_est, scale = sigma_est, shape = xi_est)
  }, error = function(e) NULL)
  
  if (!is.null(lmom_start)) {
    fit <- tryCatch(
      evd::fgev(bmx, start = lmom_start),
      error = function(e) NULL
    )
    if (.valid(fit)) return(fit)
  }
  
  # ------------------------------------------------------------------
  # Attempt 3: Conservative Gumbel-like start (ξ ≈ 0)
  # ------------------------------------------------------------------
  gumbel_start <- tryCatch({
    sigma_g <- sd(bmx) * sqrt(6) / pi
    mu_g    <- mean(bmx) - 0.5772 * sigma_g
    list(loc = mu_g, scale = sigma_g, shape = 0.01)
  }, error = function(e) NULL)
  
  if (!is.null(gumbel_start)) {
    fit <- tryCatch(
      evd::fgev(bmx, start = gumbel_start),
      error = function(e) NULL
    )
    if (.valid(fit)) return(fit)
  }
  
  # ------------------------------------------------------------------
  # Attempt 4: Shape-clamped Nelder-Mead (derivative-free)
  #   BFGS fails when the likelihood surface is flat or ridged (common
  #   with near-constant block maxima from collinear/sparse DGPs).
  #   Nelder-Mead is more tolerant of these geometries. We clamp ξ to
  #   [-0.5, 0.5] via a penalty term to avoid degenerate Fréchet/Weibull
  #   tails that make the likelihood unbounded.
  # ------------------------------------------------------------------
  
  # Penalised negative log-likelihood with shape clamp
  # (defined here so both Attempt 4 and 5 can use it)
  nll_gev_clamped <- function(par, data) {
    mu    <- par[1]
    sigma <- par[2]
    xi    <- par[3]
    
    # Hard constraints: sigma > 0, |xi| <= 0.5
    if (sigma <= 1e-10) return(1e12)
    if (abs(xi) > 0.5)  return(1e12)
    
    n <- length(data)
    
    if (abs(xi) < 1e-8) {
      # Gumbel case (ξ → 0)
      z <- (data - mu) / sigma
      nll <- n * log(sigma) + sum(z) + sum(exp(-z))
    } else {
      z <- 1 + xi * (data - mu) / sigma
      # All z must be > 0 for the GEV density to be defined
      if (any(z <= 0)) return(1e12)
      nll <- n * log(sigma) + (1 + 1/xi) * sum(log(z)) + sum(z^(-1/xi))
    }
    
    if (!is.finite(nll)) return(1e12)
    return(nll)
  }
  
  nm_fit <- tryCatch({
    
    # Use the best available starting values (Gumbel or L-moments)
    if (!is.null(gumbel_start)) {
      init <- c(gumbel_start$loc, gumbel_start$scale, gumbel_start$shape)
    } else if (!is.null(lmom_start)) {
      init <- c(lmom_start$loc, lmom_start$scale, lmom_start$shape)
    } else {
      init <- c(median(bmx), max(sd(bmx), 1e-4), 0.01)
    }
    
    optres <- optim(
      par = init,
      fn  = nll_gev_clamped,
      data = bmx,
      method  = "Nelder-Mead",
      control = list(maxit = 5000, reltol = 1e-8)
    )
    
    if (optres$convergence != 0) stop("Nelder-Mead did not converge")
    if (optres$par[2] <= 0)      stop("Negative scale from Nelder-Mead")
    
    # Wrap into an fgev-compatible list so downstream code works unchanged
    list(
      estimate = c(loc = optres$par[1], scale = optres$par[2], shape = optres$par[3]),
      deviance = 2 * optres$value,
      convergence = optres$convergence,
      data = bmx
    )
  }, error = function(e) NULL)
  
  if (!is.null(nm_fit) && nm_fit$estimate["scale"] > 0) return(nm_fit)
  
  # ------------------------------------------------------------------
  # Attempt 5: Median/IQR robust start with Nelder-Mead
  #   When the block maxima are heavily contaminated or near-constant,
  #   mean/sd-based initialisations place the optimizer in a dead zone.
  #   Median and IQR are breakdown-resistant (50% breakdown point) and
  #   give a reasonable starting region even for degenerate samples.
  # ------------------------------------------------------------------
  robust_nm_fit <- tryCatch({
    med_bmx <- median(bmx)
    iqr_bmx <- max(IQR(bmx), 1e-6)  # floor to prevent zero scale
    
    # IQR ≈ 1.573σ for Gumbel, so σ ≈ IQR / 1.573
    sigma_r <- iqr_bmx / 1.573
    mu_r    <- med_bmx - 0.3665 * sigma_r  # Gumbel median ≈ μ - σ·ln(ln2)
    
    init_r <- c(mu_r, sigma_r, 0.0)
    
    optres_r <- optim(
      par = init_r,
      fn  = nll_gev_clamped,
      data = bmx,
      method  = "Nelder-Mead",
      control = list(maxit = 5000, reltol = 1e-8)
    )
    
    if (optres_r$convergence != 0) stop("Robust NM did not converge")
    if (optres_r$par[2] <= 0)      stop("Negative scale from robust NM")
    
    list(
      estimate = c(loc = optres_r$par[1], scale = optres_r$par[2], shape = optres_r$par[3]),
      deviance = 2 * optres_r$value,
      convergence = optres_r$convergence,
      data = bmx
    )
  }, error = function(e) NULL)
  
  if (!is.null(robust_nm_fit) && robust_nm_fit$estimate["scale"] > 0) {
    return(robust_nm_fit)
  }
  
  # All attempts failed
  return(NULL)
}