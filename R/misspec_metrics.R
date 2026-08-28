# ==============================================================================
# File: R/misspec_metrics.R
# Purpose: Evaluation metrics and detection-boundary summaries for Script 05.
# ==============================================================================

compute_residual_abnormality05 <- function(
    model,
    cutoff = 2.5
) {
  
  r <-
    stats::residuals(
      model
    )
  
  
  r_med <-
    stats::median(
      r,
      na.rm = TRUE
    )
  
  
  # Robust residual scale.
  #
  # stats::mad(..., constant = 1.4826)
  # already applies the normal-consistency correction.
  residual_mad <-
    stats::mad(
      r,
      center = r_med,
      constant = 1.4826,
      na.rm = TRUE
    )
  
  
  # Safe fallback
  if (
    !is.finite(residual_mad) ||
    residual_mad <= .Machine$double.eps
  ) {
    
    residual_mad <-
      stats::sd(
        r,
        na.rm = TRUE
      )
  }
  
  
  if (
    !is.finite(residual_mad) ||
    residual_mad <= .Machine$double.eps
  ) {
    
    return(
      list(
        
        # Existing-style outputs
        outlier_rate =
          NA_real_,
        
        n_abnormal =
          NA_integer_,
        
        indices =
          integer(0),
        
        # Compatibility aliases
        rate =
          NA_real_,
        
        n =
          NA_integer_,
        
        # New outputs
        residual_mad =
          NA_real_,
        
        q95_abs_std_resid =
          NA_real_,
        
        max_abs_std_resid =
          NA_real_
      )
    )
  }
  
  
  std_resid <-
    (r - r_med) /
    residual_mad
  
  
  abs_std_resid <-
    abs(
      std_resid
    )
  
  
  abnormal_idx <-
    which(
      abs_std_resid >
        cutoff
    )
  
  
  outlier_rate <-
    length(abnormal_idx) /
    length(r)
  
  
  n_abnormal <-
    length(
      abnormal_idx
    )
  
  
  q95_abs_std_resid <-
    as.numeric(
      stats::quantile(
        abs_std_resid,
        probs = 0.95,
        na.rm = TRUE,
        names = FALSE
      )
    )
  
  
  max_abs_std_resid <-
    max(
      abs_std_resid,
      na.rm = TRUE
    )
  
  
  list(
    
    # Existing outputs
    outlier_rate =
      outlier_rate,
    
    n_abnormal =
      n_abnormal,
    
    indices =
      abnormal_idx,
    
    # Compatibility aliases in case old engine uses $rate / $n
    rate =
      outlier_rate,
    
    n =
      n_abnormal,
    
    # New outputs
    residual_mad =
      residual_mad,
    
    q95_abs_std_resid =
      q95_abs_std_resid,
    
    max_abs_std_resid =
      max_abs_std_resid
  )
}

compute_selection_overlap05 <- function(selected_idx, affected_idx) {
  if (is.null(affected_idx) || length(affected_idx) == 0L) return(NA_real_)
  if (length(selected_idx) == 0L) return(0)
  length(intersect(selected_idx, affected_idx)) / length(affected_idx)
}

compute_precision05 <- function(selected_idx, affected_idx) {
  if (is.null(affected_idx) || length(affected_idx) == 0L) return(NA_real_)
  if (length(selected_idx) == 0L) return(NA_real_)
  length(intersect(selected_idx, affected_idx)) / length(selected_idx)
}

compute_estimation_metrics05 <- function(
    coef,
    se,
    true_beta
) {
  coef <-
    as.numeric(coef)[1L]
  
  se <-
    as.numeric(se)[1L]
  
  true_beta <-
    as.numeric(true_beta)[1L]
  
  
  bias <- if (
    is.finite(coef) &&
    is.finite(true_beta)
  ) {
    
    coef - true_beta
    
  } else {
    
    NA_real_
  }
  
  
  abs_bias <-
    if (is.finite(bias)) {
      abs(bias)
    } else {
      NA_real_
    }
  
  
  coverage <- if (
    is.finite(coef) &&
    is.finite(se) &&
    se >= 0 &&
    is.finite(true_beta)
  ) {
    
    as.integer(
      true_beta >= coef - 1.96 * se &&
        true_beta <= coef + 1.96 * se
    )
    
  } else {
    
    NA_integer_
  }
  
  
  c(
    bias = bias,
    abs_bias = abs_bias,
    coverage = coverage
  )
}

# Post-processing helper: derive a null cutoff from severity_level == 0 and then
# detection probability at each positive severity level.
compute_detection_power05 <- function(results, statistic = "standardized_delta",
                                      alpha = 0.05,
                                      group_vars = c("n", "x_type", "error_type",
                                                     "scenario", "k_fraction", "diagnostic")) {
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stat_sym <- rlang::sym(statistic)
  
  nulls <- results |>
    dplyr::filter(.data$severity_level == 0L, .data$model_state == "wrong") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(null_cutoff = stats::quantile(!!stat_sym, 1 - alpha, na.rm = TRUE), .groups = "drop")
  
  results |>
    dplyr::filter(.data$model_state == "wrong") |>
    dplyr::left_join(nulls, by = group_vars) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_vars, "severity_level", "severity_target")))) |>
    dplyr::summarise(
      detection_power = mean((!!stat_sym) > .data$null_cutoff, na.rm = TRUE),
      .groups = "drop"
    )
}

estimate_detection_boundary05 <- function(power_table, target_power = 0.80,
                                          group_vars = c("n", "x_type", "error_type",
                                                         "scenario", "k_fraction", "diagnostic")) {
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  power_table |>
    dplyr::filter(.data$severity_level > 0L, .data$detection_power >= target_power) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::slice_min(.data$severity_target, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::rename(detection_boundary = severity_target)
}