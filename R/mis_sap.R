# ==============================================================================
# MIS-SAP: selection-adjusted permutation influential-coalition analysis
#
# Core contract:
#   1. Search both influence directions over a bounded multiscale k-grid.
#   2. Repeat the complete search under every residual permutation.
#   3. Use the global permutation p-value for coalition detection.
#   4. After rejection, exactly refit the target coefficient after deleting the
#      calibrated sensitivity coalition.
#   5. Optionally select between the two directional coalitions using a robust
#      coefficient anchor.
#
# The robust anchor is strictly post-detection and cannot alter the statistic,
# p-value, rejection decision, selected k, or sensitivity coalition.
#
# No global rejection, invalid anchor, invalid refit, or directional tie implies
# no automatic deletion recommendation.
# ==============================================================================


#' MIS-SAP selection-adjusted influential-coalition analysis
#'
#' @param formula Linear-model formula with one finite numeric response.
#' @param data Data frame containing all variables required by `formula`.
#' @param target Optional exact model-matrix coefficient name.
#' @param target_pos Optional integer model-matrix coefficient position.
#'   When both target arguments are supplied, they must identify the same
#'   non-intercept coefficient.
#' @param k_grid Optional positive integer coalition-size grid. `NULL` uses the
#'   bounded adaptive grid.
#' @param B_perm Positive integer number of residual permutations.
#' @param alpha Global significance level in `(0, 1)`.
#' @param max_fraction Maximum candidate coalition fraction.
#' @param use_robust_anchor Whether to compute the optional post-rejection
#'   robust-cleaning recommendation.
#' @param anchor_coefficient Optional finite robust target-coefficient estimate.
#'   `NULL` uses `robustbase::lmrob()`.
#' @param tie_tolerance Finite nonnegative numerical comparison tolerance. It is
#'   not a statistical uncertainty radius.
#'
#' @return Object of class `mis_sap`.
#' @export


# ------------------------------------------------------------------------------
# 0. Strict internal validators
# ------------------------------------------------------------------------------

.sap_assert_flag <- function(x, name) {
  if (
    !is.logical(x) ||
    length(x) != 1L ||
    is.na(x)
  ) {
    stop(
      name,
      " must be one non-missing logical value.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


.sap_assert_number <- function(
    x,
    name,
    lower = -Inf,
    upper = Inf,
    lower_open = FALSE,
    upper_open = FALSE
) {
  scalar_numeric <- (
    is.numeric(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      is.finite(x)
  )
  
  valid_lower <- scalar_numeric && if (lower_open) {
    x > lower
  } else {
    x >= lower
  }
  
  valid_upper <- scalar_numeric && if (upper_open) {
    x < upper
  } else {
    x <= upper
  }
  
  if (
    !scalar_numeric ||
    !valid_lower ||
    !valid_upper
  ) {
    left_bracket <- if (lower_open) "(" else "["
    right_bracket <- if (upper_open) ")" else "]"
    
    stop(
      name,
      " must be one finite numeric value in ",
      left_bracket,
      lower,
      ", ",
      upper,
      right_bracket,
      ".",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


.sap_assert_positive_integer <- function(x, name) {
  if (
    !is.numeric(x) ||
    length(x) != 1L ||
    is.na(x) ||
    !is.finite(x) ||
    x != as.integer(x) ||
    x < 1L
  ) {
    stop(
      name,
      " must be one positive integer.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


.sap_validate_test_controls <- function(
    B_perm,
    alpha,
    max_fraction,
    use_robust_anchor,
    tie_tolerance
) {
  .sap_assert_positive_integer(
    B_perm,
    "`B_perm`"
  )
  
  .sap_assert_number(
    alpha,
    "`alpha`",
    lower = 0,
    upper = 1,
    lower_open = TRUE,
    upper_open = TRUE
  )
  
  .sap_assert_number(
    max_fraction,
    "`max_fraction`",
    lower = 0,
    upper = 1,
    lower_open = TRUE
  )
  
  .sap_assert_flag(
    use_robust_anchor,
    "`use_robust_anchor`"
  )
  
  .sap_assert_number(
    tie_tolerance,
    "`tie_tolerance`",
    lower = 0
  )
  
  B_perm <- as.integer(
    B_perm
  )
  
  minimum_attainable_p <- 1 / (
    B_perm +
      1L
  )
  
  if (minimum_attainable_p > alpha) {
    stop(
      "`B_perm` is too small for rejection at the requested `alpha`: ",
      "the minimum attainable p-value is ",
      format(
        minimum_attainable_p,
        digits = 8
      ),
      ".",
      call. = FALSE
    )
  }
  
  invisible(
    list(
      B_perm = B_perm,
      alpha = alpha,
      max_fraction = max_fraction,
      use_robust_anchor = use_robust_anchor,
      tie_tolerance = tie_tolerance,
      minimum_attainable_p = minimum_attainable_p
    )
  )
}


.sap_validate_delete_indices <- function(
    delete,
    n
) {
  if (is.null(delete)) {
    return(integer(0L))
  }
  
  if (length(delete) == 0L) {
    return(integer(0L))
  }
  
  if (
    is.factor(delete) ||
    is.character(delete) ||
    is.logical(delete) ||
    !is.numeric(delete) ||
    anyNA(delete) ||
    any(!is.finite(delete)) ||
    any(delete != as.integer(delete)) ||
    any(delete < 1L) ||
    any(delete > n)
  ) {
    stop(
      "`delete` must contain valid integer row positions.",
      call. = FALSE
    )
  }
  
  sort(
    unique(
      as.integer(delete)
    )
  )
}


.sap_validate_target_position <- function(
    target_pos,
    coefficient_names
) {
  if (
    !is.numeric(target_pos) ||
    length(target_pos) != 1L ||
    is.na(target_pos) ||
    !is.finite(target_pos) ||
    target_pos != as.integer(target_pos) ||
    target_pos < 1L ||
    target_pos > length(coefficient_names)
  ) {
    stop(
      "`target_pos` is outside the model-matrix column range.",
      call. = FALSE
    )
  }
  
  as.integer(target_pos)
}


# ------------------------------------------------------------------------------
# 1. Model construction and target resolution
# ------------------------------------------------------------------------------

.sap_build_model_components <- function(
    formula,
    data
) {
  if (!inherits(formula, "formula")) {
    stop(
      "`formula` must be a model formula.",
      call. = FALSE
    )
  }
  
  if (
    !is.data.frame(data) ||
    nrow(data) < 3L
  ) {
    stop(
      "`data` must be a data frame with at least three rows.",
      call. = FALSE
    )
  }
  
  model_frame <- stats::model.frame(
    formula = formula,
    data = data,
    na.action = stats::na.fail,
    drop.unused.levels = TRUE
  )
  
  terms_object <- stats::terms(
    model_frame
  )
  
  if (length(attr(terms_object, "offset")) > 0L) {
    stop(
      "MIS-SAP currently does not support model offsets.",
      call. = FALSE
    )
  }
  
  response <- stats::model.response(
    model_frame
  )
  
  if (
    !is.numeric(response) ||
    is.matrix(response) ||
    anyNA(response) ||
    any(!is.finite(response))
  ) {
    stop(
      "MIS-SAP currently requires one finite numeric response.",
      call. = FALSE
    )
  }
  
  model_matrix <- stats::model.matrix(
    terms_object,
    data = model_frame
  )
  
  if (
    !is.numeric(model_matrix) ||
    anyNA(model_matrix) ||
    any(!is.finite(model_matrix))
  ) {
    stop(
      "The model matrix must contain only finite numeric values.",
      call. = FALSE
    )
  }
  
  if (
    nrow(model_matrix) <=
    ncol(model_matrix) +
    1L
  ) {
    stop(
      "The fitted design leaves insufficient residual degrees of freedom.",
      call. = FALSE
    )
  }
  
  coefficient_names <- colnames(
    model_matrix
  )
  
  if (
    is.null(coefficient_names) ||
    anyNA(coefficient_names) ||
    any(!nzchar(coefficient_names))
  ) {
    coefficient_names <- paste0(
      "coefficient_",
      seq_len(
        ncol(model_matrix)
      )
    )
  }
  
  row_positions <- match(
    row.names(model_frame),
    row.names(data)
  )
  
  if (
    anyNA(row_positions) ||
    length(unique(row_positions)) !=
    nrow(model_frame)
  ) {
    row_positions <- seq_len(
      nrow(model_frame)
    )
  }
  
  list(
    formula = formula,
    data = data,
    model_frame = model_frame,
    terms = terms_object,
    X = model_matrix,
    y = as.numeric(response),
    coefficient_names = coefficient_names,
    model_row_names = row.names(model_frame),
    data_row_positions = as.integer(row_positions),
    n = nrow(model_matrix),
    p = ncol(model_matrix)
  )
}


.sap_resolve_target <- function(
    components,
    target,
    target_pos
) {
  coefficient_names <- components$coefficient_names
  
  target_supplied <- !is.null(target)
  position_supplied <- !is.null(target_pos)
  
  if (!target_supplied && !position_supplied) {
    stop(
      "Supply either `target` or `target_pos`.",
      call. = FALSE
    )
  }
  
  if (target_supplied) {
    if (
      !is.character(target) ||
      length(target) != 1L ||
      is.na(target) ||
      !nzchar(target)
    ) {
      stop(
        "`target` must be one nonempty coefficient name.",
        call. = FALSE
      )
    }
    
    matches <- which(
      coefficient_names ==
        target
    )
    
    if (length(matches) != 1L) {
      stop(
        "`target` must exactly identify one model-matrix coefficient. ",
        "Available names are: ",
        paste(
          coefficient_names,
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }
    
    resolved_position <- as.integer(
      matches
    )
  } else {
    resolved_position <- .sap_validate_target_position(
      target_pos = target_pos,
      coefficient_names = coefficient_names
    )
  }
  
  if (target_supplied && position_supplied) {
    validated_position <- .sap_validate_target_position(
      target_pos = target_pos,
      coefficient_names = coefficient_names
    )
    
    if (!identical(resolved_position, validated_position)) {
      stop(
        "`target` and `target_pos` identify different coefficients.",
        call. = FALSE
      )
    }
  }
  
  if (
    identical(
      coefficient_names[
        resolved_position
      ],
      "(Intercept)"
    )
  ) {
    stop(
      "The intercept cannot be used as the MIS-SAP target coefficient.",
      call. = FALSE
    )
  }
  
  list(
    name = coefficient_names[
      resolved_position
    ],
    position = resolved_position
  )
}


# ------------------------------------------------------------------------------
# 2. FWL projection and exact deleted refitting
# ------------------------------------------------------------------------------

.sap_fwl_components <- function(
    X,
    y,
    target_pos,
    energy_tolerance = sqrt(.Machine$double.eps)
) {
  target_pos <- .sap_validate_target_position(
    target_pos = target_pos,
    coefficient_names = colnames(X)
  )
  
  target_predictor <- X[
    ,
    target_pos
  ]
  
  if (ncol(X) == 1L) {
    y_fwl <- as.numeric(y)
    x_fwl <- as.numeric(target_predictor)
    nuisance_rank <- 0L
  } else {
    nuisance <- X[
      ,
      -target_pos,
      drop = FALSE
    ]
    
    nuisance_qr <- qr(
      nuisance,
      LAPACK = FALSE
    )
    
    if (nuisance_qr$rank < ncol(nuisance)) {
      stop(
        "The nuisance design is rank deficient.",
        call. = FALSE
      )
    }
    
    if (exists("fwl", mode = "function", inherits = TRUE)) {
      projected <- fwl(
        y = y,
        X = target_predictor,
        Z = nuisance
      )
      
      if (
        !is.matrix(projected) ||
        ncol(projected) < 2L
      ) {
        stop(
          "`fwl()` did not return the expected two-column projection.",
          call. = FALSE
        )
      }
      
      y_fwl <- as.numeric(
        projected[
          ,
          1L
        ]
      )
      
      x_fwl <- as.numeric(
        projected[
          ,
          2L
        ]
      )
    } else {
      y_fwl <- as.numeric(
        qr.resid(
          nuisance_qr,
          y
        )
      )
      
      x_fwl <- as.numeric(
        qr.resid(
          nuisance_qr,
          target_predictor
        )
      )
    }
    
    nuisance_rank <- nuisance_qr$rank
  }
  
  target_energy <- sum(
    x_fwl^2
  )
  
  scale_reference <- max(
    1,
    sum(
      target_predictor^2
    )
  )
  
  if (
    !is.finite(target_energy) ||
    target_energy <=
    energy_tolerance *
    scale_reference
  ) {
    stop(
      "The retained target predictor has insufficient residualized energy.",
      call. = FALSE
    )
  }
  
  coefficient <- sum(
    x_fwl *
      y_fwl
  ) /
    target_energy
  
  fitted <- x_fwl *
    coefficient
  
  residual <- y_fwl -
    fitted
  
  list(
    y_fwl = y_fwl,
    x_fwl = x_fwl,
    coefficient = unname(coefficient),
    fitted = fitted,
    residual = residual,
    target_energy = unname(target_energy),
    nuisance_rank = as.integer(nuisance_rank)
  )
}


.sap_exact_refit_target <- function(
    components,
    target_pos,
    delete = integer(0L),
    energy_tolerance = sqrt(.Machine$double.eps)
) {
  target_pos <- .sap_validate_target_position(
    target_pos = target_pos,
    coefficient_names = components$coefficient_names
  )
  
  delete <- .sap_validate_delete_indices(
    delete = delete,
    n = components$n
  )
  
  keep <- setdiff(
    seq_len(
      components$n
    ),
    delete
  )
  
  if (
    length(keep) <=
    components$p
  ) {
    stop(
      "The deletion leaves insufficient observations for refitting.",
      call. = FALSE
    )
  }
  
  full_fit <- .sap_fwl_components(
    X = components$X,
    y = components$y,
    target_pos = target_pos,
    energy_tolerance = energy_tolerance
  )
  
  retained_fit <- .sap_fwl_components(
    X = components$X[
      keep,
      ,
      drop = FALSE
    ],
    y = components$y[
      keep
    ],
    target_pos = target_pos,
    energy_tolerance = energy_tolerance
  )
  
  denominator_fraction <-
    retained_fit$target_energy /
    full_fit$target_energy
  
  if (
    !is.finite(denominator_fraction) ||
    denominator_fraction <= 0
  ) {
    stop(
      "The retained exact-refit denominator is not positive.",
      call. = FALSE
    )
  }
  
  list(
    exact = TRUE,
    target_name =
      components$coefficient_names[
        target_pos
      ],
    target_pos = target_pos,
    deleted_set = delete,
    retained_set = as.integer(keep),
    deleted_n = length(delete),
    retained_n = length(keep),
    coefficient_full =
      full_fit$coefficient,
    coefficient_retained =
      retained_fit$coefficient,
    exact_shift =
      retained_fit$coefficient -
      full_fit$coefficient,
    denominator_full =
      full_fit$target_energy,
    denominator_retained =
      retained_fit$target_energy,
    denominator_fraction =
      unname(
        denominator_fraction
      )
  )
}


# ------------------------------------------------------------------------------
# 3. Bounded adaptive k-grid
# ------------------------------------------------------------------------------

.sap_adaptive_k_grid <- function(
    n,
    design_columns,
    fixed_small = c(
      1L,
      2L,
      5L,
      10L,
      20L
    ),
    proportions = c(
      0.005,
      0.010,
      0.025,
      0.050
    ),
    max_fraction = 0.05
) {
  .sap_assert_positive_integer(
    n,
    "`n`"
  )
  
  .sap_assert_positive_integer(
    design_columns,
    "`design_columns`"
  )
  
  .sap_assert_number(
    max_fraction,
    "`max_fraction`",
    lower = 0,
    upper = 1,
    lower_open = TRUE
  )
  
  if (
    !is.numeric(fixed_small) ||
    anyNA(fixed_small) ||
    any(!is.finite(fixed_small)) ||
    any(fixed_small != as.integer(fixed_small)) ||
    any(fixed_small < 1L)
  ) {
    stop(
      "`fixed_small` must contain positive finite integers.",
      call. = FALSE
    )
  }
  
  if (
    !is.numeric(proportions) ||
    length(proportions) < 1L ||
    anyNA(proportions) ||
    any(!is.finite(proportions)) ||
    any(proportions <= 0) ||
    any(proportions > max_fraction)
  ) {
    stop(
      "`proportions` must contain finite values in (0, max_fraction].",
      call. = FALSE
    )
  }
  
  fraction_cap <- floor(
    n *
      max_fraction
  )
  
  degrees_of_freedom_cap <-
    n -
    design_columns -
    1L
  
  maximum_admissible <- min(
    fraction_cap,
    degrees_of_freedom_cap
  )
  
  if (maximum_admissible < 1L) {
    stop(
      "No positive candidate coalition size is admissible.",
      call. = FALSE
    )
  }
  
  proportional_candidates <- floor(
    n *
      proportions
  )
  
  grid <- sort(
    unique(
      as.integer(
        c(
          fixed_small,
          proportional_candidates
        )
      )
    )
  )
  
  grid <- grid[
    grid >= 1L &
      grid <= maximum_admissible
  ]
  
  if (length(grid) == 0L) {
    stop(
      "No adaptive candidate coalition size remains after validation.",
      call. = FALSE
    )
  }
  
  grid
}


.sap_validate_explicit_k_grid <- function(
    k_grid,
    n,
    design_columns,
    max_fraction
) {
  if (
    !is.numeric(k_grid) ||
    length(k_grid) < 1L ||
    anyNA(k_grid) ||
    any(!is.finite(k_grid)) ||
    any(k_grid != as.integer(k_grid))
  ) {
    stop(
      "`k_grid` must contain finite integer-valued candidates.",
      call. = FALSE
    )
  }
  
  fraction_cap <- floor(
    n *
      max_fraction
  )
  
  degrees_of_freedom_cap <-
    n -
    design_columns -
    1L
  
  maximum_admissible <- min(
    fraction_cap,
    degrees_of_freedom_cap
  )
  
  grid <- sort(
    unique(
      as.integer(
        k_grid
      )
    )
  )
  
  grid <- grid[
    grid >= 1L &
      grid <= maximum_admissible
  ]
  
  if (length(grid) == 0L) {
    stop(
      "No explicit candidate size remains after applying the ",
      "fraction and degrees-of-freedom limits.",
      call. = FALSE
    )
  }
  
  grid
}


# ------------------------------------------------------------------------------
# 4. Multiscale Dinkelbach search
# ------------------------------------------------------------------------------

.sap_require_search_dependencies <- function() {
  if (
    !exists(
      "dinkelbach_topk",
      mode = "function",
      inherits = TRUE
    )
  ) {
    stop(
      "`dinkelbach_topk()` is unavailable. Source `R/dinkelbach_topk.R` ",
      "or load the package namespace first.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

.sap_profile_dinkelbach <- function(
    x, r, k_grid, sgn, sum_x2,
    max_iter = 50L, tol = 1e-9
) {
  n_val <- sgn * (x * r)
  d_val <- -(x^2)
  K     <- length(k_grid)
  
  lambda    <- rep(0, K)
  idx_list  <- vector("list", K)
  converged <- rep(FALSE, K)
  
  eval_profile <- function(lambda_ref) {
    w     <- n_val - lambda_ref * d_val
    ord   <- order(w, decreasing = TRUE)
    cum_n <- cumsum(n_val[ord])[k_grid]
    cum_d <- cumsum(d_val[ord])[k_grid]
    list(ratio = cum_n / (sum_x2 + cum_d), ord = ord)
  }
  
  for (j in seq_len(K)) {
    if (j > 1L && !converged[j]) {
      lambda[j] <- max(lambda[j], lambda[j - 1L])
    }
    ev <- NULL
    for (iter in seq_len(max_iter)) {
      ev    <- eval_profile(lambda[j])
      new_l <- ev$ratio[j]
      hint <- !converged & ev$ratio > lambda
      hint[j] <- FALSE
      lambda[hint] <- ev$ratio[hint]
      if (abs(new_l - lambda[j]) < tol) {
        lambda[j] <- new_l
        break
      }
      lambda[j] <- new_l
    }
    idx_list[[j]] <- ev$ord[seq_len(k_grid[j])]
    converged[j]  <- TRUE
  }
  
  list(lambda = lambda, dfbeta = sgn * lambda, indices = idx_list)
}

.sap_search_multiscale <- function(
    x,
    r,
    k_grid
) {
  .sap_require_search_dependencies()
  
  if (
    !is.numeric(x) ||
    !is.numeric(r) ||
    length(x) != length(r) ||
    length(x) < 3L ||
    anyNA(x) ||
    anyNA(r) ||
    any(!is.finite(x)) ||
    any(!is.finite(r))
  ) {
    stop(
      "`x` and `r` must be finite numeric vectors of equal length.",
      call. = FALSE
    )
  }
  
  sum_x2 <- sum(
    x^2
  )
  
  if (
    !is.finite(sum_x2) ||
    sum_x2 <= 0
  ) {
    stop(
      "The residualized target-predictor energy must be positive.",
      call. = FALSE
    )
  }
  
  prof_pos <- .sap_profile_dinkelbach(
    x = x, r = r, k_grid = k_grid, sgn =  1L, sum_x2 = sum_x2
  )
  prof_neg <- .sap_profile_dinkelbach(
    x = x, r = r, k_grid = k_grid, sgn = -1L, sum_x2 = sum_x2
  )
  
  searches <- lapply(
    seq_along(k_grid),
    function(j) {
      shift_pos <- prof_pos$dfbeta[j]
      shift_neg <- prof_neg$dfbeta[j]
      if (abs(shift_pos) >= abs(shift_neg)) {
        selected_direction <- 1L
        selected_statistic <- abs(shift_pos)
        selected_indices   <- prof_pos$indices[[j]]
      } else {
        selected_direction <- -1L
        selected_statistic <- abs(shift_neg)
        selected_indices   <- prof_neg$indices[[j]]
      }
      list(
        k         = as.integer(k_grid[j]),
        direction = as.integer(selected_direction),
        statistic = unname(selected_statistic),
        indices   = as.integer(selected_indices),
        positive  = list(
          indices = as.integer(prof_pos$indices[[j]]),
          shift   = unname(shift_pos)
        ),
        negative  = list(
          indices = as.integer(prof_neg$indices[[j]]),
          shift   = unname(shift_neg)
        )
      )
    }
  )
  
  list(
    profile = data.frame(
      k = vapply(
        searches,
        `[[`,
        integer(1),
        "k"
      ),
      direction = vapply(
        searches,
        `[[`,
        integer(1),
        "direction"
      ),
      observed_stat = vapply(
        searches,
        `[[`,
        numeric(1),
        "statistic"
      ),
      stringsAsFactors = FALSE
    ),
    sets = lapply(
      searches,
      `[[`,
      "indices"
    ),
    searches = searches
  )
}


.sap_directional_candidates_at_k <- function(
    fwl_fit,
    k
) {
  search <- .sap_search_multiscale(
    x = fwl_fit$x_fwl,
    r = fwl_fit$residual,
    k_grid = as.integer(k)
  )
  
  fit <- search$searches[[
    1L
  ]]
  
  list(
    positive_set =
      fit$positive$indices,
    negative_set =
      fit$negative$indices,
    positive_search_shift =
      fit$positive$shift,
    negative_search_shift =
      fit$negative$shift,
    sensitivity_direction =
      fit$direction,
    sensitivity_set =
      fit$indices
  )
}


# ------------------------------------------------------------------------------
# 5. Selection-adjusted permutation calibration
# ------------------------------------------------------------------------------

.sap_calibrated_test <- function(
    components,
    target_pos,
    k_grid,
    B_perm,
    alpha,
    max_fraction
) {
  if (is.null(k_grid)) {
    effective_grid <- .sap_adaptive_k_grid(
      n = components$n,
      design_columns = components$p,
      max_fraction = max_fraction
    )
    
    grid_source <- "bounded_adaptive"
  } else {
    effective_grid <- .sap_validate_explicit_k_grid(
      k_grid = k_grid,
      n = components$n,
      design_columns = components$p,
      max_fraction = max_fraction
    )
    
    grid_source <- "explicit"
  }
  
  fwl_fit <- .sap_fwl_components(
    X = components$X,
    y = components$y,
    target_pos = target_pos
  )
  
  observed <- .sap_search_multiscale(
    x = fwl_fit$x_fwl,
    r = fwl_fit$residual,
    k_grid = effective_grid
  )
  
  perm_stats <- matrix(
    NA_real_,
    nrow = B_perm,
    ncol = length(effective_grid)
  )
  
  observed_global_now <- max(observed$profile$observed_stat)
  
  # Besag–Clifford early stop: once `h` exceedances of the global max
  # are seen, p >= h/(b+1) > alpha is guaranteed -> cannot reject.
  h_stop      <- ceiling(alpha * (B_perm + 1))
  exceedances <- 0L
  b_used      <- 0L
  
  for (b in seq_len(B_perm)) {
    y_permuted <- fwl_fit$fitted +
      sample(fwl_fit$residual, replace = FALSE)
    
    beta_permuted <- sum(fwl_fit$x_fwl * y_permuted) /
      fwl_fit$target_energy
    
    residual_permuted <- y_permuted -
      fwl_fit$x_fwl * beta_permuted
    
    permuted_search <- .sap_search_multiscale(
      x      = fwl_fit$x_fwl,
      r      = residual_permuted,
      k_grid = effective_grid
    )
    
    perm_stats[b, ] <- permuted_search$profile$observed_stat
    b_used <- b
    
    if (max(perm_stats[b, ]) >= observed_global_now) {
      exceedances <- exceedances + 1L
      if (exceedances >= h_stop) break
    }
  }
  
  perm_stats <- perm_stats[seq_len(b_used), , drop = FALSE]
  
  null_mean <- colMeans(
    perm_stats
  )
  
  null_q95 <- apply(
    perm_stats,
    2L,
    stats::quantile,
    probs = 0.95,
    names = FALSE
  )
  
  p_by_k <- vapply(
    seq_along(
      effective_grid
    ),
    function(j) {
      (
        1 +
          sum(
            perm_stats[
              ,
              j
            ] >=
              observed$profile$observed_stat[
                j
              ]
          )
      ) /
        (
          b_used +
            1
        )
    },
    numeric(1)
  )
  
  excess_ratio <-
    observed$profile$observed_stat /
    null_q95
  
  observed_global <- max(
    observed$profile$observed_stat
  )
  
  perm_global <- apply(
    perm_stats,
    1L,
    max
  )
  
  global_p <- (
    1 +
      sum(
        perm_global >=
          observed_global
      )
  ) /
    (
      b_used +
        1
    )
  
  profile <- observed$profile
  profile$null_mean <- null_mean
  profile$null_q95 <- null_q95
  profile$p_selection_adjusted <- p_by_k
  profile$excess_ratio <- excess_ratio
  
  if (global_p <= alpha) {
    best_scale <- which.max(
      excess_ratio
    )
    
    detected <- TRUE
    selected_k <- effective_grid[
      best_scale
    ]
    selected_direction <- profile$direction[
      best_scale
    ]
    selected_set <- observed$sets[[
      best_scale
    ]]
    peak_excess_ratio <- excess_ratio[
      best_scale
    ]
    stop_reason <- "significant"
  } else {
    detected <- FALSE
    selected_k <- 0L
    selected_direction <- 0L
    selected_set <- integer(0L)
    peak_excess_ratio <- max(
      excess_ratio
    )
    stop_reason <- "not_significant"
  }
  
  list(
    detected = detected,
    selected_set = as.integer(
      selected_set
    ),
    selected_k = as.integer(
      selected_k
    ),
    selected_direction = as.integer(
      selected_direction
    ),
    global_p = unname(
      global_p
    ),
    observed_global_stat = unname(
      observed_global
    ),
    null_global_q95 = unname(
      stats::quantile(
        perm_global,
        0.95,
        names = FALSE
      )
    ),
    peak_excess_ratio = unname(
      peak_excess_ratio
    ),
    profile = profile,
    permutation_count = as.integer(
      b_used
    ),
    min_attainable_p = 1 / (
      b_used +
        1
    ),
    effective_k_grid = as.integer(
      effective_grid
    ),
    grid_source = grid_source,
    fwl_fit = fwl_fit,
    observed_search = observed,
    converged = TRUE,
    stop_reason = stop_reason
  )
}


# ------------------------------------------------------------------------------
# 6. Post-rejection robust-anchor direction
# ------------------------------------------------------------------------------

.sap_choose_direction_from_anchor <- function(
    positive_coefficient,
    negative_coefficient,
    anchor_coefficient,
    tie_tolerance = 1e-10
) {
  finite_scalar <- function(x) {
    is.numeric(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      is.finite(x)
  }
  
  .sap_assert_number(
    tie_tolerance,
    "`tie_tolerance`",
    lower = 0
  )
  
  anchor_valid <- finite_scalar(
    anchor_coefficient
  )
  
  candidate_valid <- finite_scalar(
    positive_coefficient
  ) &&
    finite_scalar(
      negative_coefficient
    )
  
  if (!anchor_valid) {
    return(
      list(
        selected_direction = 0L,
        positive_distance = NA_real_,
        negative_distance = NA_real_,
        anchor_valid = FALSE,
        candidate_valid = candidate_valid,
        tied = FALSE,
        decision_status =
          "invalid_anchor_no_deletion"
      )
    )
  }
  
  if (!candidate_valid) {
    return(
      list(
        selected_direction = 0L,
        positive_distance = NA_real_,
        negative_distance = NA_real_,
        anchor_valid = TRUE,
        candidate_valid = FALSE,
        tied = FALSE,
        decision_status =
          "invalid_candidate_no_deletion"
      )
    )
  }
  
  positive_distance <- abs(
    positive_coefficient -
      anchor_coefficient
  )
  
  negative_distance <- abs(
    negative_coefficient -
      anchor_coefficient
  )
  
  comparison_scale <- max(
    1,
    abs(
      positive_coefficient
    ),
    abs(
      negative_coefficient
    ),
    abs(
      anchor_coefficient
    )
  )
  
  tied <- abs(
    positive_distance -
      negative_distance
  ) <=
    tie_tolerance *
    comparison_scale
  
  if (tied) {
    selected_direction <- 0L
    decision_status <- "anchor_tie_no_deletion"
  } else if (
    positive_distance <
    negative_distance
  ) {
    selected_direction <- 1L
    decision_status <- "positive_direction"
  } else {
    selected_direction <- -1L
    decision_status <- "negative_direction"
  }
  
  list(
    selected_direction = as.integer(
      selected_direction
    ),
    positive_distance = unname(
      positive_distance
    ),
    negative_distance = unname(
      negative_distance
    ),
    anchor_valid = TRUE,
    candidate_valid = TRUE,
    tied = tied,
    decision_status = decision_status
  )
}


.sap_no_cleaning_result <- function(
    status,
    anchor_source
) {
  list(
    selected_direction = 0L,
    selected_set = integer(0L),
    selected_coefficient = NA_real_,
    positive_set = integer(0L),
    negative_set = integer(0L),
    positive_coefficient = NA_real_,
    negative_coefficient = NA_real_,
    anchor_coefficient = NA_real_,
    positive_distance = NA_real_,
    negative_distance = NA_real_,
    anchor_source = anchor_source,
    anchor_valid = FALSE,
    candidate_valid = FALSE,
    tied = FALSE,
    automatic_deletion_permitted = FALSE,
    decision_status = status,
    anchor_error = NA_character_,
    positive_error = NA_character_,
    negative_error = NA_character_
  )
}


.sap_robust_anchor_direction <- function(
    components,
    target_pos,
    positive_set,
    negative_set,
    anchor_coefficient = NULL,
    tie_tolerance = 1e-10
) {
  positive_set <- .sap_validate_delete_indices(
    delete = positive_set,
    n = components$n
  )
  
  negative_set <- .sap_validate_delete_indices(
    delete = negative_set,
    n = components$n
  )
  
  if (
    length(positive_set) == 0L ||
    length(negative_set) == 0L
  ) {
    stop(
      "Both directional candidates must contain at least one row.",
      call. = FALSE
    )
  }
  
  if (
    length(positive_set) !=
    length(negative_set)
  ) {
    stop(
      "The two directional candidates must have the same size.",
      call. = FALSE
    )
  }
  
  anchor_source <- if (
    is.null(anchor_coefficient)
  ) {
    "robustbase::lmrob"
  } else {
    "supplied"
  }
  
  anchor_error <- NA_character_
  
  if (is.null(anchor_coefficient)) {
    if (
      !requireNamespace(
        "robustbase",
        quietly = TRUE
      )
    ) {
      result <- .sap_no_cleaning_result(
        status =
          "invalid_anchor_no_deletion",
        anchor_source =
          anchor_source
      )
      
      result$anchor_error <-
        "Package `robustbase` is unavailable."
      
      result$positive_set <- positive_set
      result$negative_set <- negative_set
      
      return(result)
    }
    
    robust_fit <- tryCatch(
      suppressWarnings(
        robustbase::lmrob(
          formula = components$formula,
          data = components$data,
          na.action = stats::na.fail
        )
      ),
      error = function(e) {
        anchor_error <<-
          conditionMessage(e)
        
        NULL
      }
    )
    
    if (is.null(robust_fit)) {
      result <- .sap_no_cleaning_result(
        status =
          "invalid_anchor_no_deletion",
        anchor_source =
          anchor_source
      )
      
      result$anchor_error <- anchor_error
      result$positive_set <- positive_set
      result$negative_set <- negative_set
      
      return(result)
    }
    
    robust_coefficients <- stats::coef(
      robust_fit
    )
    
    target_name <-
      components$coefficient_names[
        target_pos
      ]
    
    if (
      !is.null(
        names(
          robust_coefficients
        )
      ) &&
      target_name %in%
      names(
        robust_coefficients
      )
    ) {
      anchor_coefficient <- unname(
        robust_coefficients[
          target_name
        ]
      )
    } else if (
      target_pos <=
      length(
        robust_coefficients
      )
    ) {
      anchor_coefficient <- unname(
        robust_coefficients[
          target_pos
        ]
      )
    } else {
      anchor_coefficient <- NA_real_
      anchor_error <-
        "The target coefficient is unavailable in the robust fit."
    }
  }
  
  anchor_valid <-
    is.numeric(anchor_coefficient) &&
    length(anchor_coefficient) == 1L &&
    !is.na(anchor_coefficient) &&
    is.finite(anchor_coefficient)
  
  if (!anchor_valid) {
    result <- .sap_no_cleaning_result(
      status =
        "invalid_anchor_no_deletion",
      anchor_source =
        anchor_source
    )
    
    result$anchor_error <- if (
      is.na(anchor_error)
    ) {
      "The robust anchor is not a finite scalar."
    } else {
      anchor_error
    }
    
    result$positive_set <- positive_set
    result$negative_set <- negative_set
    
    return(result)
  }
  
  positive_error <- NA_character_
  negative_error <- NA_character_
  
  positive_result <- tryCatch(
    .sap_exact_refit_target(
      components = components,
      target_pos = target_pos,
      delete = positive_set
    ),
    error = function(e) {
      positive_error <<-
        conditionMessage(e)
      
      NULL
    }
  )
  
  negative_result <- tryCatch(
    .sap_exact_refit_target(
      components = components,
      target_pos = target_pos,
      delete = negative_set
    ),
    error = function(e) {
      negative_error <<-
        conditionMessage(e)
      
      NULL
    }
  )
  
  if (
    is.null(positive_result) ||
    is.null(negative_result)
  ) {
    result <- .sap_no_cleaning_result(
      status =
        "invalid_candidate_no_deletion",
      anchor_source =
        anchor_source
    )
    
    result$anchor_coefficient <- unname(
      anchor_coefficient
    )
    result$anchor_valid <- TRUE
    result$positive_set <- positive_set
    result$negative_set <- negative_set
    result$positive_error <- positive_error
    result$negative_error <- negative_error
    
    return(result)
  }
  
  decision <- .sap_choose_direction_from_anchor(
    positive_coefficient =
      positive_result$
      coefficient_retained,
    negative_coefficient =
      negative_result$
      coefficient_retained,
    anchor_coefficient =
      anchor_coefficient,
    tie_tolerance =
      tie_tolerance
  )
  
  if (
    decision$selected_direction ==
    1L
  ) {
    selected_set <- positive_set
    selected_coefficient <-
      positive_result$
      coefficient_retained
  } else if (
    decision$selected_direction ==
    -1L
  ) {
    selected_set <- negative_set
    selected_coefficient <-
      negative_result$
      coefficient_retained
  } else {
    selected_set <- integer(0L)
    selected_coefficient <- NA_real_
  }
  
  list(
    selected_direction =
      decision$selected_direction,
    selected_set =
      as.integer(
        selected_set
      ),
    selected_coefficient =
      unname(
        selected_coefficient
      ),
    positive_set =
      positive_set,
    negative_set =
      negative_set,
    positive_coefficient =
      positive_result$
      coefficient_retained,
    negative_coefficient =
      negative_result$
      coefficient_retained,
    anchor_coefficient =
      unname(
        anchor_coefficient
      ),
    positive_distance =
      decision$positive_distance,
    negative_distance =
      decision$negative_distance,
    positive_denominator_fraction =
      positive_result$
      denominator_fraction,
    negative_denominator_fraction =
      negative_result$
      denominator_fraction,
    anchor_source =
      anchor_source,
    anchor_valid =
      decision$anchor_valid,
    candidate_valid =
      decision$candidate_valid,
    tied =
      decision$tied,
    automatic_deletion_permitted =
      decision$selected_direction !=
      0L,
    decision_status =
      decision$decision_status,
    anchor_error =
      anchor_error,
    positive_error =
      positive_error,
    negative_error =
      negative_error
  )
}


# ------------------------------------------------------------------------------
# 7. Public MIS-SAP interface
# ------------------------------------------------------------------------------

#' MIS-SAP influential-coalition analysis
#'
#' Runs a selection-adjusted permutation test over a bounded multiscale
#' coalition-size grid. If the global test rejects, the function reports the
#' calibrated sensitivity coalition and its exact deleted-refit coefficient
#' effect. It may also report an optional robust-anchor cleaning coalition.
#'
#' The robust anchor is evaluated only after global rejection. It does not alter
#' the global statistic, permutation distribution, p-value, rejection decision,
#' selected coalition size, sensitivity direction, or sensitivity coalition.
#'
#' @param formula Linear-model formula.
#' @param data Data frame containing all model variables.
#' @param target Optional exact model-matrix coefficient name, such as `"x"`.
#' @param target_pos Optional integer model-matrix coefficient position. Supply
#'   this only when `target` is not used, or supply both consistently.
#' @param k_grid Optional positive integer coalition-size grid. `NULL` uses the
#'   bounded adaptive grid.
#' @param B_perm Positive integer number of residual permutations. Use 199 for
#'   simulation development and 999 or more for final reported applications.
#' @param alpha Global significance level in `(0, 1)`.
#' @param max_fraction Maximum candidate coalition fraction. The validated
#'   prospective default is 0.05.
#' @param use_robust_anchor Whether to compute the optional post-rejection
#'   robust-cleaning direction.
#' @param anchor_coefficient Optional externally supplied finite robust anchor.
#'   `NULL` uses `robustbase::lmrob()`.
#' @param tie_tolerance Nonnegative numerical tolerance for declaring the two
#'   robust-direction candidates tied.
#'
#' @return An object of class `mis_sap`.
#' @export
mis_sap <- function(
    formula,
    data,
    target = NULL,
    target_pos = NULL,
    k_grid = NULL,
    B_perm = 999L,
    alpha = 0.05,
    max_fraction = 0.05,
    use_robust_anchor = TRUE,
    anchor_coefficient = NULL,
    tie_tolerance = 1e-10
) {
  controls <- .sap_validate_test_controls(
    B_perm = B_perm,
    alpha = alpha,
    max_fraction = max_fraction,
    use_robust_anchor = use_robust_anchor,
    tie_tolerance = tie_tolerance
  )
  
  components <- .sap_build_model_components(
    formula = formula,
    data = data
  )
  
  resolved_target <- .sap_resolve_target(
    components = components,
    target = target,
    target_pos = target_pos
  )
  
  calibrated <- .sap_calibrated_test(
    components = components,
    target_pos =
      resolved_target$position,
    k_grid = k_grid,
    B_perm = controls$B_perm,
    alpha = controls$alpha,
    max_fraction =
      controls$max_fraction
  )
  
  if (!calibrated$detected) {
    state <- "no_excessive_influence"
    
    sensitivity_exact <- NULL
    directional_candidates <- NULL
    
    cleaning <- .sap_no_cleaning_result(
      status =
        "not_detected_no_deletion",
      anchor_source =
        if (use_robust_anchor) {
          "not_evaluated"
        } else {
          "disabled"
        }
    )
  } else {
    directional_candidates <-
      .sap_directional_candidates_at_k(
        fwl_fit =
          calibrated$fwl_fit,
        k =
          calibrated$selected_k
      )
    
    reconstruction_consistent <-
      identical(
        sort(
          calibrated$selected_set
        ),
        sort(
          directional_candidates$
            sensitivity_set
        )
      ) &&
      identical(
        calibrated$
          selected_direction,
        directional_candidates$
          sensitivity_direction
      )
    
    if (!reconstruction_consistent) {
      stop(
        "The post-detection directional reconstruction does not match ",
        "the calibrated sensitivity result.",
        call. = FALSE
      )
    }
    
    sensitivity_exact <- .sap_exact_refit_target(
      components = components,
      target_pos =
        resolved_target$position,
      delete =
        calibrated$selected_set
    )
    
    if (use_robust_anchor) {
      cleaning <- .sap_robust_anchor_direction(
        components = components,
        target_pos =
          resolved_target$position,
        positive_set =
          directional_candidates$
          positive_set,
        negative_set =
          directional_candidates$
          negative_set,
        anchor_coefficient =
          anchor_coefficient,
        tie_tolerance =
          controls$tie_tolerance
      )
    } else {
      cleaning <- .sap_no_cleaning_result(
        status =
          "robust_anchor_disabled",
        anchor_source =
          "disabled"
      )
      
      cleaning$positive_set <-
        directional_candidates$
        positive_set
      
      cleaning$negative_set <-
        directional_candidates$
        negative_set
    }
    
    state <- if (
      isTRUE(
        cleaning$
        automatic_deletion_permitted
      )
    ) {
      "robust_cleaning_available"
    } else {
      "sensitivity_only"
    }
  }
  
  sensitivity_set_model <- if (
    calibrated$detected
  ) {
    calibrated$selected_set
  } else {
    integer(0L)
  }
  
  sensitivity_set_data <-
    components$data_row_positions[
      sensitivity_set_model
    ]
  
  cleaning_set_model <- cleaning$selected_set
  
  cleaning_set_data <-
    components$data_row_positions[
      cleaning_set_model
    ]
  
  result <- list(
    method = "MIS-SAP",
    version = "2.0.0",
    state = state,
    
    target = list(
      name =
        resolved_target$name,
      position =
        resolved_target$position
    ),
    
    global = list(
      reject =
        calibrated$detected,
      p_value =
        calibrated$global_p,
      alpha =
        controls$alpha,
      statistic =
        calibrated$
        observed_global_stat,
      null_q95 =
        calibrated$
        null_global_q95,
      peak_excess_ratio =
        calibrated$
        peak_excess_ratio,
      permutations =
        calibrated$
        permutation_count,
      minimum_attainable_p =
        calibrated$
        min_attainable_p,
      stop_reason =
        calibrated$
        stop_reason
    ),
    
    grid = list(
      source =
        calibrated$
        grid_source,
      values =
        calibrated$
        effective_k_grid,
      max_fraction =
        controls$
        max_fraction
    ),
    
    sensitivity = list(
      detected =
        calibrated$
        detected,
      selected_k =
        calibrated$
        selected_k,
      direction =
        calibrated$
        selected_direction,
      set =
        as.integer(
          sensitivity_set_model
        ),
      data_rows =
        as.integer(
          sensitivity_set_data
        ),
      row_names =
        components$
        model_row_names[
          sensitivity_set_model
        ],
      exact =
        sensitivity_exact
    ),
    
    cleaning = list(
      permitted =
        cleaning$
        automatic_deletion_permitted,
      direction =
        cleaning$
        selected_direction,
      set =
        as.integer(
          cleaning_set_model
        ),
      data_rows =
        as.integer(
          cleaning_set_data
        ),
      row_names =
        components$
        model_row_names[
          cleaning_set_model
        ],
      coefficient_retained =
        cleaning$
        selected_coefficient,
      anchor_coefficient =
        cleaning$
        anchor_coefficient,
      anchor_source =
        cleaning$
        anchor_source,
      anchor_valid =
        cleaning$
        anchor_valid,
      candidate_valid =
        cleaning$
        candidate_valid,
      tied =
        cleaning$
        tied,
      status =
        cleaning$
        decision_status,
      positive_set =
        cleaning$
        positive_set,
      negative_set =
        cleaning$
        negative_set,
      positive_coefficient =
        cleaning$
        positive_coefficient,
      negative_coefficient =
        cleaning$
        negative_coefficient,
      positive_distance =
        cleaning$
        positive_distance,
      negative_distance =
        cleaning$
        negative_distance,
      anchor_error =
        cleaning$
        anchor_error,
      positive_error =
        cleaning$
        positive_error,
      negative_error =
        cleaning$
        negative_error
    ),
    
    directional_candidates =
      directional_candidates,
    
    profile =
      calibrated$
      profile,
    
    assumptions = list(
      permutation = paste(
        "Residual exchangeability under the null is required for the",
        "permutation calibration."
      ),
      observations = paste(
        "Observations must be independent or use an appropriately",
        "justified exchangeability structure."
      ),
      exact_refit = paste(
        "Exact deleted refitting requires a full-rank retained nuisance",
        "design and positive retained target-predictor energy."
      ),
      cleaning = paste(
        "Robust cleaning is a post-rejection directional recommendation,",
        "not a guarantee of contamination recovery or valid",
        "post-selection confidence intervals."
      )
    ),
    
    call = match.call()
  )
  
  class(result) <- c(
    "mis_sap",
    "list"
  )
  
  result
}


# ------------------------------------------------------------------------------
# 8. S3 print method
# ------------------------------------------------------------------------------

#' @export
print.mis_sap <- function(
    x,
    ...
) {
  format_set <- function(index_set) {
    if (length(index_set) == 0L) {
      return("<empty>")
    }
    
    paste(
      index_set,
      collapse = ", "
    )
  }
  
  format_number <- function(value, digits = 5L) {
    if (
      !is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !is.finite(value)
    ) {
      return("<undefined>")
    }
    
    format(
      value,
      digits = digits
    )
  }
  
  cat("MIS-SAP influential-coalition analysis\n")
  cat("----------------------------------------\n")
  cat("State:", x$state, "\n")
  cat("Target:", x$target$name, "\n")
  cat(
    "Global p-value:",
    format_number(
      x$global$p_value,
      digits = 4L
    ),
    "\n"
  )
  cat(
    "Global rejection:",
    if (isTRUE(x$global$reject)) {
      "yes"
    } else {
      "no"
    },
    "\n"
  )
  cat(
    "Permutation count:",
    x$global$permutations,
    "\n"
  )
  cat(
    "Candidate k-grid:",
    format_set(
      x$grid$values
    ),
    "\n"
  )
  
  if (!isTRUE(x$global$reject)) {
    cat(
      "No influential coalition was detected; no deletion is recommended.\n"
    )
    
    return(
      invisible(x)
    )
  }
  
  cat(
    "Selected sensitivity k:",
    x$sensitivity$selected_k,
    "\n"
  )
  cat(
    "Sensitivity direction:",
    x$sensitivity$direction,
    "\n"
  )
  cat(
    "Sensitivity set:",
    format_set(
      x$sensitivity$set
    ),
    "\n"
  )
  cat(
    "Full-sample target coefficient:",
    format_number(
      x$sensitivity$exact$
        coefficient_full
    ),
    "\n"
  )
  cat(
    "Retained-sample target coefficient:",
    format_number(
      x$sensitivity$exact$
        coefficient_retained
    ),
    "\n"
  )
  cat(
    "Exact coefficient shift:",
    format_number(
      x$sensitivity$exact$
        exact_shift
    ),
    "\n"
  )
  
  cat(
    "Robust cleaning recommendation:",
    if (isTRUE(x$cleaning$permitted)) {
      "available"
    } else {
      "not available"
    },
    "\n"
  )
  
  if (isTRUE(x$cleaning$permitted)) {
    cat(
      "Robust-cleaning direction:",
      x$cleaning$direction,
      "\n"
    )
    cat(
      "Robust-cleaning set:",
      format_set(
        x$cleaning$set
      ),
      "\n"
    )
    cat(
      "Robust anchor coefficient:",
      format_number(
        x$cleaning$anchor_coefficient
      ),
      "\n"
    )
    cat(
      "Coefficient after robust cleaning:",
      format_number(
        x$cleaning$coefficient_retained
      ),
      "\n"
    )
  } else {
    cat(
      "Cleaning status:",
      x$cleaning$status,
      "\n"
    )
  }
  
  invisible(x)
}