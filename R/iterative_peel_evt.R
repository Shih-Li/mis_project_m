# ==============================================================================
# File: /R/iterative_peel_evt.R
# Purpose: EVT-guided iterative peel for MIS detection. Replaces the three
#          ad-hoc heuristics in iterative_peel_v2 (leverage_k gate, sigma
#          direction, sigma_tol stopping) with a single principled criterion:
#          the EVT p-value from Konrad & Kuschnig (2025).
#
# Dependencies: dinkelbach_topk.R, evt_peel.R (replaces evt_iter_dm),
#               exact_dfb_bmx.R, helpers_local.R,
#               estimators_robust.R (fit_clean_ols)
# ==============================================================================


#' EVT-Guided Iterative Peel MIS Detection
#'
#' Removes influential observations iteratively, using the EVT p-value
#' as the sole criterion for direction selection, outlier detection,
#' and stopping. Each removal is justified by a formal test of whether
#' the current most-influential set's effect is excessive under the
#' GEV null distribution.
#'
#' @param formula    A formula (e.g., y ~ x).
#' @param data       A data.frame containing the variables.
#' @param target_var Character; name of the target coefficient (e.g., "x").
#'        Used for coefficient extraction in trajectory recording.
#' @param target_pos Integer; position of the target coefficient in the
#'        design matrix (default = 2, first slope with intercept).
#' @param batch_size Integer; number of points to peel per iteration.
#'        1 = most careful (tests each point individually), higher values
#'        trade precision for speed. Default = 1.
#' @param max_iter   Integer; hard ceiling on peel iterations. Default = 100.
#'        With EVT stopping, this should rarely bind.
#' @param max_k_frac Numeric; hard cap on total fraction of n removed.
#'        Default = 0.10 (10%). Safety valve to prevent catastrophic
#'        over-removal if the EVT test is miscalibrated.
#' @param alpha      Numeric; significance level for the EVT test.
#'        Default = 0.05. At each iteration, if the most significant
#'        direction has p >= alpha, peeling stops.
#' @param block_count Integer; requested number of blocks for EVT block
#'        maxima. Adaptively reduced if the clean sample is too small.
#'        Default = 20.
#' @param min_obs_per_block Integer; minimum observations per block to
#'        ensure stable Dinkelbach solutions within each block.
#'        Default = 30.
#' @param verbose    Logical; print iteration trace. Default = FALSE.
#'
#' @return A list with components:
#'   \item{excluded}{Integer vector of all removed observation indices
#'                   (in original data row numbering).}
#'   \item{k_total}{Integer; total number removed.}
#'   \item{n_iters}{Integer; number of peel iterations performed.}
#'   \item{stop_reason}{Character; why peeling stopped. One of:
#'         "not_significant", "evt_failed", "max_k_reached",
#'         "max_iter", "insufficient_blocks", "insufficient_data".}
#'   \item{beta_trajectory}{Numeric vector; OLS coefficient after each
#'         peel step (on full data minus all excluded).}
#'   \item{sigma_trajectory}{Numeric vector; residual SE after each step
#'         (diagnostic only — not used for stopping).}
#'   \item{pval_trajectory}{Numeric vector; EVT p-value of the removed
#'         set at each step. Final element may be >= alpha (the
#'         non-significant test that triggered stopping).}
#'   \item{direction_trajectory}{Integer vector; +1 or -1 indicating
#'         which influence direction was more significant at each step.}
#'   \item{converged_trajectory}{Logical vector; whether the EVT fit
#'         converged at each step.}
#' @export
iterative_peel_evt <- function(formula, data,
                               target_var  = "x",
                               target_pos  = 2L,
                               batch_size  = 1L,
                               max_iter    = 100L,
                               max_k_frac  = 0.10,
                               alpha       = 0.05,
                               block_count = 20L,
                               min_obs_per_block = 30L,
                               verbose     = FALSE) {
  
  n_total   <- nrow(data)
  max_k_abs <- floor(n_total * max_k_frac)
  
  # ------------------------------------------------------------------
  # Response variable name (for extracting y from clean data)
  # ------------------------------------------------------------------
  response_var <- all.vars(formula[[2]])[1]
  
  # ------------------------------------------------------------------
  # State tracking
  # ------------------------------------------------------------------
  excluded_all      <- integer(0)
  beta_traj         <- numeric(0)
  sigma_traj        <- numeric(0)
  pval_traj         <- numeric(0)
  direction_traj    <- integer(0)
  converged_traj    <- logical(0)
  stop_reason       <- "max_iter"
  it                <- 0L
  
  for (it in seq_len(max_iter)) {
    
    # ----------------------------------------------------------------
    # Guard: hard cap on total removals
    # ----------------------------------------------------------------
    if (length(excluded_all) >= max_k_abs) {
      stop_reason <- "max_k_reached"
      break
    }
    
    # ----------------------------------------------------------------
    # A. Current clean data
    # ----------------------------------------------------------------
    if (length(excluded_all) > 0) {
      clean_data <- data[-excluded_all, , drop = FALSE]
    } else {
      clean_data <- data
    }
    
    n_clean <- nrow(clean_data)
    
    # Minimum sample check: need enough for OLS + meaningful EVT
    X_check <- tryCatch(
      stats::model.matrix(formula, data = clean_data),
      error = function(e) NULL
    )
    if (is.null(X_check) || n_clean < ncol(X_check) + 10L) {
      stop_reason <- "insufficient_data"
      break
    }
    
    # ----------------------------------------------------------------
    # B. Adaptive block count for current clean sample
    #
    #    Constraint: each block must hold >= batch_size observations
    #    for the Dinkelbach solver, AND >= min_obs_per_block for
    #    stable block maxima. GEV needs >= 3 blocks.
    # ----------------------------------------------------------------
    n_for_blocks <- n_clean - batch_size  # observations available after set removal
    B_max <- floor(n_for_blocks / max(batch_size, min_obs_per_block))
    B_eff <- min(block_count, B_max)
    
    if (B_eff < 3L) {
      stop_reason <- "insufficient_blocks"
      break
    }
    
    # ----------------------------------------------------------------
    # C. Fit OLS on clean data → find candidates in both directions
    # ----------------------------------------------------------------
    mod_ols_clean <- stats::lm(formula, data = clean_data)
    
    this_batch <- min(batch_size, max_k_abs - length(excluded_all))
    if (this_batch < 1L) {
      stop_reason <- "max_k_reached"
      break
    }
    
    idx_pos <- dinkelbach_topk_lm(mod_ols_clean, pos = target_pos,
                                  sign =  1L, k = this_batch)
    idx_neg <- dinkelbach_topk_lm(mod_ols_clean, pos = target_pos,
                                  sign = -1L, k = this_batch)
    
    # ----------------------------------------------------------------
    # D. Extract FWL components from clean data for EVT
    #
    #    x_target = the target predictor column (pos)
    #    Z        = all other columns (intercept + controls)
    # ----------------------------------------------------------------
    X_clean    <- stats::model.matrix(formula, data = clean_data)
    y_clean    <- clean_data[[response_var]]
    x_target   <- X_clean[, target_pos]
    Z_clean    <- X_clean[, -target_pos, drop = FALSE]
    
    # ----------------------------------------------------------------
    # E. EVT test in both directions
    #
    #    Each call runs the full pipeline: FWL → DFBETA → exact block
    #    maxima → marginal tail analysis → constrained GEV + robust
    #    fallback → M-scaling → p-value.
    #    Uses evt_peel which combines all four components.
    # ----------------------------------------------------------------
    run_evt <- function(set_indices) {
      tryCatch({
        evt_peel(
          y = y_clean, x = x_target, Z = Z_clean,
          set = set_indices, block_count = B_eff
        )
      }, error = function(e) {
        data.frame(
          shape = NA_real_, scale = NA_real_, loc = NA_real_,
          set_dfb = NA_real_, p_value = NA_real_, converged = FALSE,
          tail_coef = NA_real_, tail_source = NA_character_,
          stringsAsFactors = FALSE
        )
      })
    }
    
    evt_pos <- run_evt(idx_pos)
    evt_neg <- run_evt(idx_neg)
    
    p_pos <- if (isTRUE(evt_pos$converged)) evt_pos$p_value else NA_real_
    p_neg <- if (isTRUE(evt_neg$converged)) evt_neg$p_value else NA_real_
    
    # ----------------------------------------------------------------
    # F. Direction selection: pick the more significant direction
    #
    #    Smaller p-value = more excessive influence = correct direction
    #    to peel. If both fail, we cannot assess significance.
    # ----------------------------------------------------------------
    if (is.na(p_pos) && is.na(p_neg)) {
      stop_reason <- "evt_failed"
      converged_traj <- c(converged_traj, FALSE)
      break
    }
    
    if (is.na(p_pos)) {
      best_p <- p_neg; best_idx <- idx_neg; best_dir <- -1L
      best_converged <- TRUE
    } else if (is.na(p_neg)) {
      best_p <- p_pos; best_idx <- idx_pos; best_dir <- 1L
      best_converged <- TRUE
    } else if (p_pos <= p_neg) {
      best_p <- p_pos; best_idx <- idx_pos; best_dir <- 1L
      best_converged <- TRUE
    } else {
      best_p <- p_neg; best_idx <- idx_neg; best_dir <- -1L
      best_converged <- TRUE
    }
    
    # ----------------------------------------------------------------
    # G. THE STOPPING CRITERION: EVT significance test
    #
    #    If the most influential set in the best direction is NOT
    #    statistically excessive (p >= alpha), we've removed all
    #    genuine outliers — the remaining influence is within
    #    natural sampling variation. Stop.
    # ----------------------------------------------------------------
    if (best_p >= alpha) {
      stop_reason <- "not_significant"
      # Record the non-significant p-value for diagnostics
      pval_traj      <- c(pval_traj, best_p)
      direction_traj <- c(direction_traj, best_dir)
      converged_traj <- c(converged_traj, best_converged)
      break
    }
    
    # ----------------------------------------------------------------
    # H. Accept this peel step: influence is statistically excessive
    # ----------------------------------------------------------------
    
    # Map clean-data indices → original-data indices
    if (length(excluded_all) > 0) {
      original_rows <- seq_len(n_total)[-excluded_all]
    } else {
      original_rows <- seq_len(n_total)
    }
    orig_idx <- original_rows[best_idx]
    
    excluded_all <- c(excluded_all, orig_idx)
    
    # Record trajectories
    res_clean <- fit_clean_ols(formula, data = data,
                               exclude_idx = excluded_all)
    beta_new  <- unname(res_clean["coef"])
    
    mod_after <- tryCatch(
      stats::lm(formula, data = data[-excluded_all, , drop = FALSE]),
      error = function(e) NULL
    )
    sigma_new <- if (!is.null(mod_after)) {
      summary(mod_after)$sigma
    } else {
      NA_real_
    }
    
    beta_traj      <- c(beta_traj, beta_new)
    sigma_traj     <- c(sigma_traj, sigma_new)
    pval_traj      <- c(pval_traj, best_p)
    direction_traj <- c(direction_traj, best_dir)
    converged_traj <- c(converged_traj, best_converged)
    
    if (verbose) {
      cat(sprintf(
        "  Peel %03d: removed %d obs (dir=%+d), k_total=%3d, "
        , it, this_batch, best_dir, length(excluded_all)
      ))
      cat(sprintf(
        "beta=%.5f, sigma=%.5f, p=%.2e\n",
        beta_new, sigma_new, best_p
      ))
    }
  }
  
  return(list(
    excluded              = excluded_all,
    k_total               = length(excluded_all),
    n_iters               = min(it, max_iter),
    stop_reason           = stop_reason,
    beta_trajectory       = beta_traj,
    sigma_trajectory      = sigma_traj,
    pval_trajectory       = pval_traj,
    direction_trajectory  = direction_traj,
    converged_trajectory  = converged_traj
  ))
}


#' Convenience Wrapper: EVT Peel → Clean OLS Result
#'
#' Runs \code{iterative_peel_evt} and returns the clean OLS coefficient + SE
#' in the same format as \code{fit_clean_ols}, for direct comparison in
#' simulations alongside MIS-alpha, MIS-oracle, MM, LTS, etc.
#'
#' @inheritParams iterative_peel_evt
#'
#' @return A named numeric vector \code{c(coef = ..., se = ...)} matching
#'         the output format of \code{fit_clean_ols}.
#' @export
peel_evt_and_fit <- function(formula, data,
                             target_var  = "x",
                             target_pos  = 2L,
                             batch_size  = 1L,
                             max_iter    = 100L,
                             max_k_frac  = 0.10,
                             alpha       = 0.05,
                             block_count = 20L,
                             min_obs_per_block = 30L) {
  
  peel_result <- iterative_peel_evt(
    formula           = formula,
    data              = data,
    target_var        = target_var,
    target_pos        = target_pos,
    batch_size        = batch_size,
    max_iter          = max_iter,
    max_k_frac        = max_k_frac,
    alpha             = alpha,
    block_count       = block_count,
    min_obs_per_block = min_obs_per_block
  )
  
  fit_clean_ols(formula, data = data, exclude_idx = peel_result$excluded)
}


#' Full Peel Result Extractor for Enriched Simulation Output
#'
#' Runs \code{iterative_peel_evt} and returns both the clean OLS result
#' AND the peel metadata (k, stop_reason, trajectory summaries) as a
#' flat named list for easy column-binding in simulation data.frames.
#'
#' @inheritParams iterative_peel_evt
#' @param true_outliers Integer vector or NULL; true injected indices for
#'        computing detection overlap. NULL if no injection (clean data).
#'
#' @return A named list with components:
#'   \item{coef}{Numeric; clean OLS coefficient.}
#'   \item{se}{Numeric; clean OLS standard error.}
#'   \item{k_peel_evt}{Integer; total observations removed.}
#'   \item{peel_evt_stop}{Character; stopping reason.}
#'   \item{peel_evt_iters}{Integer; number of iterations.}
#'   \item{peel_evt_overlap}{Numeric; fraction of true outliers found
#'         (NA if true_outliers is NULL).}
#'   \item{peel_evt_final_p}{Numeric; p-value at the final iteration.}
#'   \item{peel_evt_sigma_ratio}{Numeric; final_sigma / initial_sigma.
#'         Values < 1 indicate the peel tightened the fit.}
#'   \item{peel_evt_min_p}{Numeric; minimum p-value across all iterations
#'         (most significant removal).}
#' @export
peel_evt_full <- function(formula, data,
                          target_var  = "x",
                          target_pos  = 2L,
                          batch_size  = 1L,
                          max_iter    = 100L,
                          max_k_frac  = 0.10,
                          alpha       = 0.05,
                          block_count = 20L,
                          min_obs_per_block = 30L,
                          true_outliers = NULL) {
  
  peel_result <- iterative_peel_evt(
    formula           = formula,
    data              = data,
    target_var        = target_var,
    target_pos        = target_pos,
    batch_size        = batch_size,
    max_iter          = max_iter,
    max_k_frac        = max_k_frac,
    alpha             = alpha,
    block_count       = block_count,
    min_obs_per_block = min_obs_per_block
  )
  
  # Clean OLS on the peeled data
  ols_result <- fit_clean_ols(formula, data = data,
                              exclude_idx = peel_result$excluded)
  
  # Detection overlap
  if (!is.null(true_outliers) && length(true_outliers) > 0) {
    overlap <- length(intersect(peel_result$excluded, true_outliers)) /
      length(true_outliers)
  } else {
    overlap <- NA_real_
  }
  
  # Trajectory summaries
  final_p <- if (length(peel_result$pval_trajectory) > 0) {
    tail(peel_result$pval_trajectory, 1)
  } else {
    NA_real_
  }
  
  min_p <- if (length(peel_result$pval_trajectory) > 0) {
    min(peel_result$pval_trajectory, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  # Sigma ratio: how much the peel tightened the fit
  if (length(peel_result$sigma_trajectory) > 0) {
    # Initial sigma (before any peeling)
    mod_init <- stats::lm(formula, data = data)
    sigma_init <- summary(mod_init)$sigma
    sigma_final <- tail(peel_result$sigma_trajectory, 1)
    sigma_ratio <- sigma_final / sigma_init
  } else {
    sigma_ratio <- 1.0  # no peeling occurred
  }
  
  list(
    coef                = unname(ols_result["coef"]),
    se                  = unname(ols_result["se"]),
    k_peel_evt          = peel_result$k_total,
    peel_evt_stop       = peel_result$stop_reason,
    peel_evt_iters      = peel_result$n_iters,
    peel_evt_overlap    = overlap,
    peel_evt_final_p    = final_p,
    peel_evt_sigma_ratio = sigma_ratio,
    peel_evt_min_p      = min_p
  )
}