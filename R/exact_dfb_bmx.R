# ==============================================================================
# File: /R/exact_dfb_bmx.R
# Purpose: Exact implementation of the block maxima DFBETA search using 
#          Dinkelbach's method for linear-fractional programming.
# Dependencies: Requires helpers_local.R (provides make_blocks)
# ==============================================================================
#'
#' Exact Set Influence Block Maxima
#'
#' @param X Vector of predictor values (full sample)
#' @param R Vector of residual values (full sample)
#' @param set Vector of indices for the influential observation set
#' @param block_count Number of blocks to divide the data into
#' @return Numeric vector of mathematically exact block maxima DFBETA values
exact_dfb_bmx <- function(X, R, set, block_count) {
  sgn <- sign(sum(X[set] * R[set]))
  if (sgn == 0) stop("dfbeta of set is exactly zero")
  
  X_inf <- X[-set]
  R_inf <- R[-set]
  
  sumX2 <- sum(X_inf^2) 
  
  nS <- length(set)
  block_size <- length(X_inf) %/% block_count
  
  Xbl <- make_blocks(X_inf, block_size)
  Rbl <- make_blocks(R_inf, block_size)
  
  res <- numeric(block_count)
  
  for (i in seq_len(block_count)) {
    x_bl <- Xbl[, i]
    r_bl <- Rbl[, i]
    
    n_val <- sgn * (x_bl * r_bl)
    d_val <- - (x_bl^2)
    
    lambda <- 0
    
    for (iter in 1:50) {
      w <- n_val - lambda * d_val
      idx <- order(w, decreasing = TRUE)[seq_len(nS)]
      new_lambda <- sum(n_val[idx]) / (sumX2 + sum(d_val[idx]))
      if (abs(new_lambda - lambda) < 1e-9) break
      lambda <- new_lambda
    }
    res[i] <- sgn * lambda
  }
  
  return(res)
}

#' Diagnostic Version of Exact Block Maxima (Dinkelbach)
#'
#' Identical computation to \code{exact_dfb_bmx}, but returns per-block
#' metadata alongside the block maxima. Used to diagnose convergence
#' failures at high k: if \code{ratio_k_to_blocksize} approaches 1.0,
#' each block's Dinkelbach solver is selecting nearly all observations,
#' producing degenerate (near-constant) block maxima that cause GEV
#' fitting to fail.
#'
#' @param X Vector of predictor values (full sample).
#' @param R Vector of residual values (full sample).
#' @param set Vector of indices for the influential observation set.
#' @param block_count Number of blocks to divide the data into.
#'
#' @return A list with components:
#'   \item{bmx}{Numeric vector of block maxima — identical to
#'              \code{exact_dfb_bmx(X, R, set, block_count)}.}
#'   \item{iterations}{Integer vector of length \code{block_count} —
#'                     Dinkelbach iterations used per block.}
#'   \item{block_size}{Integer — observations per block.}
#'   \item{effective_k}{Integer — k used within each block
#'                      (equals \code{length(set)}).}
#'   \item{ratio_k_to_blocksize}{Numeric — \code{effective_k / block_size}.
#'         Values above ~0.5 indicate degenerate block-level search.}
#' @export
exact_dfb_bmx_diag <- function(X, R, set, block_count) {
  sgn <- sign(sum(X[set] * R[set]))
  if (sgn == 0) stop("dfbeta of set is exactly zero")
  
  X_inf <- X[-set]
  R_inf <- R[-set]
  
  sumX2 <- sum(X_inf^2)
  
  nS <- length(set)
  block_size <- length(X_inf) %/% block_count
  
  Xbl <- make_blocks(X_inf, block_size)
  Rbl <- make_blocks(R_inf, block_size)
  
  res        <- numeric(block_count)
  iter_counts <- integer(block_count)
  
  for (i in seq_len(block_count)) {
    x_bl <- Xbl[, i]
    r_bl <- Rbl[, i]
    
    n_val <- sgn * (x_bl * r_bl)
    d_val <- -(x_bl^2)
    
    lambda <- 0
    
    for (iter in 1:50) {
      w <- n_val - lambda * d_val
      idx <- order(w, decreasing = TRUE)[seq_len(nS)]
      new_lambda <- sum(n_val[idx]) / (sumX2 + sum(d_val[idx]))
      if (abs(new_lambda - lambda) < 1e-9) break
      lambda <- new_lambda
    }
    
    res[i]         <- sgn * lambda
    iter_counts[i] <- iter
  }
  
  list(
    bmx                 = res,
    iterations          = iter_counts,
    block_size          = as.integer(block_size),
    effective_k         = as.integer(nS),
    ratio_k_to_blocksize = nS / block_size
  )
}