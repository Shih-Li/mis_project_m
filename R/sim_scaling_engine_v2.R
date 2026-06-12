# ==============================================================================
# File: /R/sim_scaling_engine_v2.R
# Purpose: Worker function for the revised Script 03 algorithmic comparison.
#          Runs the 2x2 detection factorial ({greedy, exact} x {no refinement,
#          with refinement}) and 3 EVT configurations ({greedy bmx + constrained,
#          exact bmx + constrained, exact bmx + robust}) per Monte Carlo draw.
#
# Dependencies (must be sourced before this file):
#   helpers_local.R, 03_scaling_dgp.R, fast_sens_topk.R, dinkelbach_topk.R,
#   exact_dfb_bmx.R, evt_iter.R, evt_iter_dm.R, metrics_influence.R
# ==============================================================================


#' Safe Outlier Injection for Scaling DGPs
#'
#' Copied from scripts/03_alg_comp.R to make the engine self-contained.
#' Shifts main covariates and twists Y to mask the target parameter.
#' @keywords internal
inject_safe_outliers_v2 <- function(dgp_res, k, magnitude = 4) {
  df <- dgp_res$data
  N  <- nrow(df)
  outlier_idx <- sample(1:N, k, replace = FALSE)
  
  shift_x <- max(mad(df$X, constant = 1.4826), 1e-4) * magnitude
  df$X[outlier_idx] <- df$X[outlier_idx] + shift_x
  
  if ("M" %in% names(df)) {
    if (length(unique(df$M)) <= 2) {
      df$M[outlier_idx] <- 1L
    } else {
      shift_m <- max(mad(df$M, constant = 1.4826), 1e-4) * magnitude
      df$M[outlier_idx] <- df$M[outlier_idx] + shift_m
    }
  }
  if ("W" %in% names(df)) {
    shift_w <- max(mad(df$W, constant = 1.4826), 1e-4) * magnitude
    df$W[outlier_idx] <- df$W[outlier_idx] + shift_w
  }
  
  z_cols <- grep("^(Z_matrix\\.|X[0-9]+|Z[0-9]+)$", names(df), value = TRUE)
  if (length(z_cols) > 0) {
    for (zc in head(z_cols, 3)) {
      shift_z <- max(mad(df[[zc]], constant = 1.4826), 1e-4) * magnitude
      df[[zc]][outlier_idx] <- df[[zc]][outlier_idx] + shift_z
    }
  }
  
  clean_mod  <- lm(dgp_res$formula, data = df[-outlier_idx, ])
  expected_y <- predict(clean_mod, newdata = df[outlier_idx, ])
  scale_y    <- max(mad(df$y, constant = 1.4826), 1e-4)
  y_shift_dir <- sample(c(1, -1), size = 1)
  df$y[outlier_idx] <- expected_y + (y_shift_dir * scale_y * magnitude * 2)
  
  dgp_res$data <- df
  dgp_res$true_outliers <- outlier_idx
  dgp_res$injection_y_direction <- y_shift_dir
  return(dgp_res)
}


#' Detect MIS with Dual-Sign Oracle Direction Selection
#'
#' Runs detection in both sign directions, picks the one with better
#' overlap against the true outlier set. Returns overlap, influence
#' ratio, and CPU time.
#'
#' @param detect_fn Function(sign) -> integer vector of indices.
#' @param true_outliers Integer vector of true injected indices.
#' @param y_fwl Numeric vector; FWL-projected response.
#' @param x_fwl Numeric vector; FWL-projected predictor.
#' @return A list with: indices, overlap, influence_ratio, cpu.
#' @keywords internal
detect_with_direction <- function(detect_fn, true_outliers, y_fwl, x_fwl, k) {
  t0 <- proc.time()[3]
  
  idx_pos <- detect_fn(1L)
  idx_neg <- detect_fn(-1L)
  
  overlap_pos <- length(intersect(idx_pos, true_outliers))
  overlap_neg <- length(intersect(idx_neg, true_outliers))
  
  if (overlap_neg > overlap_pos) {
    best_idx <- idx_neg
    best_overlap <- overlap_neg
  } else {
    best_idx <- idx_pos
    best_overlap <- overlap_pos
  }
  
  cpu <- proc.time()[3] - t0
  
  ir <- tryCatch(
    influence_ratio_fwl(y_fwl, x_fwl, best_idx, true_outliers),
    error = function(e) NA_real_
  )
  
  list(
    indices         = best_idx,
    overlap         = best_overlap / k,
    influence_ratio = ir,
    cpu             = cpu
  )
}


#' Single Iteration of the Revised Algorithmic Comparison
#'
#' One Monte Carlo draw: generates data, injects outliers, runs all 4
#' detection arms and 3 EVT configurations, returns a flat 1-row data.frame.
#'
#' @param iter_id    Integer; iteration index.
#' @param N          Integer; sample size.
#' @param k          Integer; number of injected outliers.
#' @param B          Integer; block count for EVT.
#' @param architecture Character; DGP architecture name.
#' @param magnitude  Numeric; injection severity.
#' @param iter_seed  Integer; seed for reproducibility.
#' @param rho        Numeric or NA; correlation for collinear_interaction.
#'
#' @return A 1-row data.frame.
#' @export
run_scaling_iteration_v2 <- function(iter_id, N, k, B, architecture,
                                     magnitude, iter_seed, rho = NA) {
  set.seed(iter_seed)
  
  # Failure row template
  fail_row <- data.frame(
    iter = iter_id,
    # Detection — 4 arms
    det_greedy_noref = NA_real_, det_greedy_ref = NA_real_,
    det_exact_noref  = NA_real_, det_exact_ref  = NA_real_,
    ir_greedy_noref  = NA_real_, ir_greedy_ref  = NA_real_,
    ir_exact_noref   = NA_real_, ir_exact_ref   = NA_real_,
    cpu_det_greedy_noref = NA_real_, cpu_det_greedy_ref = NA_real_,
    cpu_det_exact_noref  = NA_real_, cpu_det_exact_ref  = NA_real_,
    # EVT — 3 configurations
    p_greedy_con  = NA_real_, conv_greedy_con  = FALSE, shape_greedy_con  = NA_real_,
    p_exact_con   = NA_real_, conv_exact_con   = FALSE, shape_exact_con   = NA_real_,
    p_exact_rob   = NA_real_, conv_exact_rob   = FALSE, shape_exact_rob   = NA_real_,
    cpu_evt_greedy_con = NA_real_, cpu_evt_exact_con = NA_real_,
    cpu_evt_exact_rob  = NA_real_,
    error_msg = NA_character_,
    stringsAsFactors = FALSE
  )
  
  tryCatch({
    
    # =================================================================
    # A. Generate Data + Inject Outliers
    # =================================================================
    dgp_clean <- generate_scaling_dgp(
      N = N, architecture = architecture,
      rho = if (!is.null(rho) && !is.na(rho)) rho else 0.8
    )
    dgp_poisoned <- inject_safe_outliers_v2(dgp_clean, k = k,
                                            magnitude = magnitude)
    
    # =================================================================
    # B. Fit Base Model + Extract FWL Components
    # =================================================================
    mod  <- lm(dgp_poisoned$formula, data = dgp_poisoned$data)
    tpos <- dgp_poisoned$target_pos
    
    X_full   <- model.matrix(dgp_poisoned$formula, dgp_poisoned$data)
    x_target <- X_full[, tpos]
    Z_fwl    <- X_full[, -tpos, drop = FALSE]
    y_vec    <- dgp_poisoned$data$y
    true_out <- dgp_poisoned$true_outliers
    
    # FWL projection for influence ratio
    if (ncol(X_full) == 1L) {
      x_fwl <- X_full[, 1]
      y_fwl <- y_vec
    } else {
      qr_Z <- qr(Z_fwl)
      x_fwl <- qr.resid(qr_Z, x_target)
      y_fwl <- qr.resid(qr_Z, y_vec)
    }
    
    # =================================================================
    # C. Detection — 4 Arms
    # =================================================================
    
    # C1: Greedy, no refinement
    r_gn <- detect_with_direction(
      detect_fn = function(sgn) fast_sens_topk(mod, pos = tpos, sign = sgn,
                                               k = k, max_refine = 0),
      true_outliers = true_out, y_fwl = y_fwl, x_fwl = x_fwl, k = k
    )
    
    # C2: Greedy, with refinement
    r_gr <- detect_with_direction(
      detect_fn = function(sgn) fast_sens_topk(mod, pos = tpos, sign = sgn,
                                               k = k, max_refine = 5),
      true_outliers = true_out, y_fwl = y_fwl, x_fwl = x_fwl, k = k
    )
    
    # C3: Exact (Dinkelbach), no refinement
    r_en <- detect_with_direction(
      detect_fn = function(sgn) dinkelbach_topk_lm(mod, pos = tpos,
                                                   sign = sgn, k = k),
      true_outliers = true_out, y_fwl = y_fwl, x_fwl = x_fwl, k = k
    )
    
    # C4: Exact (Dinkelbach), with refinement
    r_er <- detect_with_direction(
      detect_fn = function(sgn) dinkelbach_topk_refined(mod, pos = tpos,
                                                        sign = sgn, k = k,
                                                        max_refine = 5),
      true_outliers = true_out, y_fwl = y_fwl, x_fwl = x_fwl, k = k
    )
    
    # =================================================================
    # D. EVT — 3 Configurations
    # =================================================================
    
    # D1: Greedy bmx + constrained GEV (on greedy_ref detected set)
    t0 <- proc.time()[3]
    evt_gc <- tryCatch(
      evt_iter(y = y_vec, x = x_target, Z = Z_fwl,
               set = r_gr$indices, block_count = B),
      error = function(e) data.frame(
        shape = NA_real_, scale = NA_real_, loc = NA_real_,
        set_dfb = NA_real_, p_value = NA_real_, converged = FALSE,
        stringsAsFactors = FALSE
      )
    )
    cpu_gc <- proc.time()[3] - t0
    
    # D2: Exact bmx + constrained GEV (on exact_ref detected set)
    t0 <- proc.time()[3]
    evt_ec <- tryCatch(
      evt_iter_dm_v2(y = y_vec, x = x_target, Z = Z_fwl,
                     set = r_er$indices, block_count = B,
                     gev_method = "constrained"),
      error = function(e) data.frame(
        shape = NA_real_, scale = NA_real_, loc = NA_real_,
        set_dfb = NA_real_, p_value = NA_real_, converged = FALSE,
        stringsAsFactors = FALSE
      )
    )
    cpu_ec <- proc.time()[3] - t0
    
    # D3: Exact bmx + robust GEV (on exact_ref detected set)
    t0 <- proc.time()[3]
    evt_er <- tryCatch(
      evt_iter_dm_v2(y = y_vec, x = x_target, Z = Z_fwl,
                     set = r_er$indices, block_count = B,
                     gev_method = "robust"),
      error = function(e) data.frame(
        shape = NA_real_, scale = NA_real_, loc = NA_real_,
        set_dfb = NA_real_, p_value = NA_real_, converged = FALSE,
        stringsAsFactors = FALSE
      )
    )
    cpu_er <- proc.time()[3] - t0
    
    # =================================================================
    # E. Assemble Output
    # =================================================================
    data.frame(
      iter = iter_id,
      # Detection overlap (fraction)
      det_greedy_noref = r_gn$overlap,
      det_greedy_ref   = r_gr$overlap,
      det_exact_noref  = r_en$overlap,
      det_exact_ref    = r_er$overlap,
      # Influence ratio
      ir_greedy_noref  = r_gn$influence_ratio,
      ir_greedy_ref    = r_gr$influence_ratio,
      ir_exact_noref   = r_en$influence_ratio,
      ir_exact_ref     = r_er$influence_ratio,
      # Detection CPU
      cpu_det_greedy_noref = r_gn$cpu,
      cpu_det_greedy_ref   = r_gr$cpu,
      cpu_det_exact_noref  = r_en$cpu,
      cpu_det_exact_ref    = r_er$cpu,
      # EVT: greedy bmx + constrained
      p_greedy_con     = evt_gc$p_value,
      conv_greedy_con  = evt_gc$converged,
      shape_greedy_con = evt_gc$shape,
      # EVT: exact bmx + constrained
      p_exact_con      = evt_ec$p_value,
      conv_exact_con   = evt_ec$converged,
      shape_exact_con  = evt_ec$shape,
      # EVT: exact bmx + robust
      p_exact_rob      = evt_er$p_value,
      conv_exact_rob   = evt_er$converged,
      shape_exact_rob  = evt_er$shape,
      # EVT CPU
      cpu_evt_greedy_con = cpu_gc,
      cpu_evt_exact_con  = cpu_ec,
      cpu_evt_exact_rob  = cpu_er,
      # Error
      error_msg = NA_character_,
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    warning(sprintf("Iter %d failed: %s", iter_id, e$message))
    fail_row$error_msg <- e$message
    fail_row
  })
}