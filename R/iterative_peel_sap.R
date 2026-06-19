# ==============================================================================
# File: /R/iterative_peel_sap.R
# Purpose: Selection-Adjusted Permutation (SAP) multiscale MIS detection.
#
#          Dinkelbach optimization deliberately selects an extreme set.
#          Therefore, a fixed-set p-value is invalid when the same data are
#          used to select and test that set. SAP corrects this post-selection
#          bias by repeating the complete Dinkelbach search under residual
#          permutation.
#
#          At each SAP test:
#            1. Residualize the target predictor and response using FWL.
#            2. Search both influence directions over a candidate k-grid.
#            3. Repeat the complete search under residual permutation.
#            4. Use the global permutation p-value for detection.
#            5. If significant, select k by the largest excess ratio:
#
#                   observed DFBETA / permutation 95% quantile.
#
#          iterative_peel_sap() optionally repeats this procedure after
#          removing the detected coalition. The validated default is one
#          detection/removal step (max_iter = 1).
#
# Dependencies:
#   - helpers_local.R   : fwl()
#   - dinkelbach_topk.R : dinkelbach_topk()
# ==============================================================================


#' Exact Multiscale Dinkelbach Search
#'
#' Runs exact Dinkelbach optimization in both coefficient directions for every
#' candidate set size and retains the direction with larger absolute DFBETA.
#'
#' @param x Numeric vector; FWL-orthogonalized target predictor.
#' @param r Numeric vector; residuals from the FWL regression.
#' @param k_grid Integer vector of candidate set sizes.
#'
#' @return A list containing a per-k profile and the selected indices at each k.
#' @keywords internal
.sap_search_multiscale <- function(x, r, k_grid) {
  
  sum_x2 <- sum(x^2)
  
  searches <- lapply(k_grid, function(k_now) {
    
    fit_pos <- dinkelbach_topk(
      x = x,
      r = r,
      k = k_now,
      sgn = 1L,
      sum_x2 = sum_x2
    )
    
    fit_neg <- dinkelbach_topk(
      x = x,
      r = r,
      k = k_now,
      sgn = -1L,
      sum_x2 = sum_x2
    )
    
    if (abs(fit_pos$dfbeta) >= abs(fit_neg$dfbeta)) {
      list(
        k = k_now,
        direction = 1L,
        statistic = abs(fit_pos$dfbeta),
        indices = fit_pos$indices
      )
    } else {
      list(
        k = k_now,
        direction = -1L,
        statistic = abs(fit_neg$dfbeta),
        indices = fit_neg$indices
      )
    }
  })
  
  list(
    profile = data.frame(
      k = vapply(searches, `[[`, integer(1), "k"),
      direction = vapply(searches, `[[`, integer(1), "direction"),
      observed_stat = vapply(searches, `[[`, numeric(1), "statistic"),
      stringsAsFactors = FALSE
    ),
    sets = lapply(searches, `[[`, "indices")
  )
}


#' Selection-Adjusted Permutation Multiscale MIS Test
#'
#' Calibrates the complete Dinkelbach search over candidate set sizes and both
#' influence directions using residual permutation.
#'
#' @param formula A model formula, for example `y ~ x`.
#' @param data A data.frame containing all model variables.
#' @param target_pos Integer; target coefficient position in the model matrix.
#' @param k_grid Integer vector of candidate influential-set sizes.
#' @param B_perm Integer; number of residual permutations.
#' @param alpha Numeric; significance level for the global test.
#' @param max_k Integer or NULL; largest candidate size allowed.
#'
#' @return A list containing the detection decision, selected set, selected k,
#'         direction, global p-value, excess ratio, and per-k profile.
#' @export
sap_multiscale_test <- function(
    formula,
    data,
    target_pos = 2L,
    k_grid = c(1L, 2L, 5L, 10L, 20L, 50L, 100L),
    B_perm = 199L,
    alpha = 0.05,
    max_k = NULL
) {
  
  mod <- stats::lm(formula, data = data)
  X <- stats::model.matrix(mod)
  y <- stats::model.response(stats::model.frame(mod))
  
  if (target_pos < 1L || target_pos > ncol(X)) {
    stop("`target_pos` is outside the model-matrix column range.")
  }
  
  if (is.null(max_k)) {
    max_k <- max(k_grid)
  }
  
  k_grid <- sort(unique(as.integer(k_grid)))
  k_grid <- k_grid[
    k_grid >= 1L &
      k_grid <= max_k &
      k_grid < nrow(X) - ncol(X)
  ]
  
  if (length(k_grid) == 0L) {
    stop("No valid candidate sizes remain in `k_grid`.")
  }
  
  # FWL projection is performed once, then reused by every search.
  if (ncol(X) == 1L) {
    Y_fwl <- y
    X_fwl <- X[, 1L]
  } else {
    fwl_vars <- fwl(
      y = y,
      X = X[, target_pos],
      Z = X[, -target_pos, drop = FALSE]
    )
    
    Y_fwl <- fwl_vars[, 1L]
    X_fwl <- fwl_vars[, 2L]
  }
  
  sum_x2 <- sum(X_fwl^2)
  beta_fwl <- sum(X_fwl * Y_fwl) / sum_x2
  fitted_fwl <- X_fwl * beta_fwl
  residual_fwl <- Y_fwl - fitted_fwl
  
  # Observed complete search.
  observed <- .sap_search_multiscale(
    x = X_fwl,
    r = residual_fwl,
    k_grid = k_grid
  )
  
  # Permutation null. Every replicate repeats the complete search.
  perm_stats <- matrix(
    NA_real_,
    nrow = B_perm,
    ncol = length(k_grid)
  )
  
  for (b in seq_len(B_perm)) {
    
    Y_perm <- fitted_fwl + sample(
      residual_fwl,
      replace = FALSE
    )
    
    beta_perm <- sum(X_fwl * Y_perm) / sum_x2
    residual_perm <- Y_perm - X_fwl * beta_perm
    
    perm_search <- .sap_search_multiscale(
      x = X_fwl,
      r = residual_perm,
      k_grid = k_grid
    )
    
    perm_stats[b, ] <- perm_search$profile$observed_stat
  }
  
  # k-specific calibration.
  null_mean <- colMeans(perm_stats)
  
  null_q95 <- apply(
    perm_stats,
    2L,
    stats::quantile,
    probs = 0.95,
    names = FALSE
  )
  
  p_by_k <- vapply(seq_along(k_grid), function(j) {
    (
      1 +
        sum(perm_stats[, j] >= observed$profile$observed_stat[j])
    ) / (
      B_perm + 1
    )
  }, numeric(1))
  
  excess_ratio <- observed$profile$observed_stat / null_q95
  
  # Global calibration repeats the search over all candidate k values.
  observed_global <- max(observed$profile$observed_stat)
  perm_global <- apply(perm_stats, 1L, max)
  
  global_p <- (
    1 +
      sum(perm_global >= observed_global)
  ) / (
    B_perm + 1
  )
  
  profile <- observed$profile
  profile$null_mean <- null_mean
  profile$null_q95 <- null_q95
  profile$p_selection_adjusted <- p_by_k
  profile$excess_ratio <- excess_ratio
  
  if (global_p <= alpha) {
    
    best_scale <- which.max(excess_ratio)
    
    detected <- TRUE
    selected_k <- k_grid[best_scale]
    selected_direction <- profile$direction[best_scale]
    selected_set <- observed$sets[[best_scale]]
    peak_excess_ratio <- excess_ratio[best_scale]
    stop_reason <- "significant"
    
  } else {
    
    detected <- FALSE
    selected_k <- 0L
    selected_direction <- 0L
    selected_set <- integer(0)
    peak_excess_ratio <- max(excess_ratio)
    stop_reason <- "not_significant"
  }
  
  list(
    detected = detected,
    selected_set = as.integer(selected_set),
    selected_k = as.integer(selected_k),
    selected_direction = as.integer(selected_direction),
    global_p = global_p,
    observed_global_stat = observed_global,
    null_global_q95 = unname(
      stats::quantile(perm_global, 0.95)
    ),
    peak_excess_ratio = peak_excess_ratio,
    profile = profile,
    permutation_count = B_perm,
    min_attainable_p = 1 / (B_perm + 1),
    converged = TRUE,
    stop_reason = stop_reason
  )
}


#' Iterative Selection-Adjusted Permutation Peel
#'
#' Applies the SAP multiscale test, removes the selected influential coalition,
#' and optionally repeats the test on the remaining observations.
#'
#' The currently validated configuration is `max_iter = 1`. When more than one
#' peel step is requested, Bonferroni alpha spending is applied across the
#' maximum number of tests.
#'
#' @param formula A model formula, for example `y ~ x`.
#' @param data A data.frame containing all model variables.
#' @param target_var Character; target coefficient name for reporting.
#' @param target_pos Integer; target coefficient position in the model matrix.
#' @param k_grid Integer vector of candidate influential-set sizes.
#' @param B_perm Integer; residual permutations per SAP test.
#' @param alpha Numeric; overall significance level.
#' @param max_iter Integer; maximum accepted peel steps. Default = 1.
#' @param max_k_frac Numeric; maximum fraction of observations removable.
#' @param verbose Logical; print iteration information.
#'
#' @return A list containing removed indices, stopping reason, trajectories,
#'         global p-values, selected sizes, directions, and excess ratios.
#' @export
iterative_peel_sap <- function(
    formula,
    data,
    target_var = "x",
    target_pos = 2L,
    k_grid = c(1L, 2L, 5L, 10L, 20L, 50L, 100L),
    B_perm = 199L,
    alpha = 0.05,
    max_iter = 1L,
    max_k_frac = 0.10,
    verbose = FALSE
) {
  
  n_total <- nrow(data)
  max_k_abs <- floor(n_total * max_k_frac)
  alpha_iter <- alpha / max_iter
  
  if (1 / (B_perm + 1) > alpha_iter) {
    stop(
      "`B_perm` is too small for the requested alpha and max_iter."
    )
  }
  
  excluded_all <- integer(0)
  beta_traj <- numeric(0)
  sigma_traj <- numeric(0)
  p_traj <- numeric(0)
  k_traj <- integer(0)
  direction_traj <- integer(0)
  excess_traj <- numeric(0)
  profile_traj <- list()
  
  stop_reason <- "max_iter"
  n_iters <- 0L
  
  for (it in seq_len(max_iter)) {
    
    remaining_cap <- max_k_abs - length(excluded_all)
    
    if (remaining_cap < 1L) {
      stop_reason <- "max_k_reached"
      break
    }
    
    original_rows <- if (length(excluded_all) == 0L) {
      seq_len(n_total)
    } else {
      seq_len(n_total)[-excluded_all]
    }
    clean_data <- data[original_rows, , drop = FALSE]
    
    sap_result <- sap_multiscale_test(
      formula = formula,
      data = clean_data,
      target_pos = target_pos,
      k_grid = k_grid,
      B_perm = B_perm,
      alpha = alpha_iter,
      max_k = remaining_cap
    )
    
    p_traj <- c(p_traj, sap_result$global_p)
    k_traj <- c(k_traj, sap_result$selected_k)
    direction_traj <- c(
      direction_traj,
      sap_result$selected_direction
    )
    excess_traj <- c(
      excess_traj,
      sap_result$peak_excess_ratio
    )
    profile_traj[[it]] <- sap_result$profile
    
    if (!sap_result$detected) {
      stop_reason <- "not_significant"
      break
    }
    
    selected_original <- original_rows[
      sap_result$selected_set
    ]
    
    excluded_all <- c(
      excluded_all,
      selected_original
    )
    
    mod_clean <- stats::lm(
      formula,
      data = data[-excluded_all, , drop = FALSE]
    )
    
    beta_new <- unname(
      stats::coef(mod_clean)[target_pos]
    )
    
    sigma_new <- summary(mod_clean)$sigma
    
    beta_traj <- c(beta_traj, beta_new)
    sigma_traj <- c(sigma_traj, sigma_new)
    n_iters <- n_iters + 1L
    
    if (verbose) {
      cat(sprintf(
        paste0(
          "  SAP %02d: removed k=%d, total=%d, ",
          "%s=%.4f, sigma=%.4f, p=%.4g, excess=%.3f\n"
        ),
        it,
        sap_result$selected_k,
        length(excluded_all),
        target_var,
        beta_new,
        sigma_new,
        sap_result$global_p,
        sap_result$peak_excess_ratio
      ))
    }
  }
  
  list(
    excluded = excluded_all,
    k_total = length(excluded_all),
    n_iters = n_iters,
    stop_reason = stop_reason,
    beta_trajectory = beta_traj,
    sigma_trajectory = sigma_traj,
    global_p_trajectory = p_traj,
    selected_k_trajectory = k_traj,
    direction_trajectory = direction_traj,
    excess_ratio_trajectory = excess_traj,
    profile_trajectory = profile_traj,
    alpha_iteration = alpha_iter,
    permutation_count = B_perm
  )
}