# ==============================================================================
# File: /R/mis_sek_adapter.R
# Purpose: Data-level adapter for independent MIS-sek certification.
#
# Public interface:
#   run_mis_sek_from_data()
#
# The adapter:
#   1. Builds the model matrix and resolves the target coefficient.
#   2. Searches both coefficient-shift directions for every admissible k.
#   3. Calibrates the two-direction profile by residual permutation.
#   4. Constructs denominator, overlap, and local-stability diagnostics.
#   5. Calls mis_sek() as an independent certification tool.
#   6. Returns a deletion set only for an automatically authorized,
#      point-identified certificate with a valid exact-refit trace.
#
# Dependencies:
#   dinkelbach_topk.R
#   mis_sek.R
#
# MIS-SAP is not called and no MIS-SAP result is consumed.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Internal validation helpers
# ------------------------------------------------------------------------------

.sek_adapter_assert_flag <- function(x, name) {
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


.sek_adapter_assert_number <- function(
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


.sek_adapter_assert_positive_integer <- function(x, name) {
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


.sek_adapter_require_dependencies <- function() {
  required <- c(
    "dinkelbach_topk",
    "mis_sek"
  )
  
  missing <- required[
    !vapply(
      required,
      exists,
      logical(1),
      mode = "function",
      inherits = TRUE
    )
  ]
  
  if (length(missing) > 0L) {
    stop(
      "Missing required function(s): ",
      paste(missing, collapse = ", "),
      ". Source `R/dinkelbach_topk.R` and `R/mis_sek.R` first.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# 1. Model construction and target resolution
# ------------------------------------------------------------------------------

.sek_adapter_build_components <- function(
    formula,
    data,
    target,
    target_pos
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
  
  terms_object <- stats::terms(model_frame)
  
  if (length(attr(terms_object, "offset")) > 0L) {
    stop(
      "MIS-sek adapter does not support model offsets.",
      call. = FALSE
    )
  }
  
  response <- stats::model.response(model_frame)
  
  if (
    !is.numeric(response) ||
    is.matrix(response) ||
    anyNA(response) ||
    any(!is.finite(response))
  ) {
    stop(
      "The response must be one finite numeric vector.",
      call. = FALSE
    )
  }
  
  X <- stats::model.matrix(
    terms_object,
    data = model_frame
  )
  
  if (
    !is.numeric(X) ||
    anyNA(X) ||
    any(!is.finite(X))
  ) {
    stop(
      "The model matrix must contain only finite numeric values.",
      call. = FALSE
    )
  }
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (n <= p + 1L) {
    stop(
      "The fitted design leaves insufficient residual degrees of freedom.",
      call. = FALSE
    )
  }
  
  coefficient_names <- colnames(X)
  
  if (
    is.null(coefficient_names) ||
    anyNA(coefficient_names) ||
    any(!nzchar(coefficient_names))
  ) {
    coefficient_names <- paste0(
      "coefficient_",
      seq_len(p)
    )
    colnames(X) <- coefficient_names
  }
  
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
    
    matches <- which(coefficient_names == target)
    
    if (length(matches) != 1L) {
      stop(
        "`target` must identify exactly one model-matrix coefficient. ",
        "Available names are: ",
        paste(coefficient_names, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    
    resolved_pos <- as.integer(matches)
  } else {
    resolved_pos <- as.integer(target_pos)
  }
  
  valid_position <- (
    is.numeric(resolved_pos) &&
      length(resolved_pos) == 1L &&
      !is.na(resolved_pos) &&
      is.finite(resolved_pos) &&
      resolved_pos == as.integer(resolved_pos) &&
      resolved_pos >= 1L &&
      resolved_pos <= p
  )
  
  if (!valid_position) {
    stop(
      "`target_pos` is outside the model-matrix column range.",
      call. = FALSE
    )
  }
  
  if (target_supplied && position_supplied) {
    supplied_pos <- as.integer(target_pos)
    
    if (
      length(supplied_pos) != 1L ||
      is.na(supplied_pos) ||
      !is.finite(supplied_pos) ||
      supplied_pos != resolved_pos
    ) {
      stop(
        "`target` and `target_pos` identify different coefficients.",
        call. = FALSE
      )
    }
  }
  
  if (
    identical(
      coefficient_names[resolved_pos],
      "(Intercept)"
    )
  ) {
    stop(
      "The intercept cannot be the MIS-sek target coefficient.",
      call. = FALSE
    )
  }
  
  row_positions <- match(
    row.names(model_frame),
    row.names(data)
  )
  
  if (
    anyNA(row_positions) ||
    length(unique(row_positions)) != n
  ) {
    row_positions <- seq_len(n)
  }
  
  full_qr <- qr(
    X,
    LAPACK = FALSE
  )
  
  if (full_qr$rank < p) {
    stop(
      "The full model matrix is rank deficient.",
      call. = FALSE
    )
  }
  
  list(
    formula = formula,
    data = data,
    model_frame = model_frame,
    terms = terms_object,
    X = X,
    y = as.numeric(response),
    n = n,
    p = p,
    coefficient_names = coefficient_names,
    target_name = coefficient_names[resolved_pos],
    target_pos = resolved_pos,
    model_row_names = row.names(model_frame),
    data_row_positions = as.integer(row_positions),
    full_qr = full_qr
  )
}


# ------------------------------------------------------------------------------
# 2. FWL projection and exact deleted refitting
# ------------------------------------------------------------------------------

.sek_adapter_fwl <- function(
    X,
    y,
    target_pos,
    energy_tolerance
) {
  target_predictor <- X[, target_pos]
  
  if (ncol(X) == 1L) {
    x_fwl <- as.numeric(target_predictor)
    y_fwl <- as.numeric(y)
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
    
    x_fwl <- as.numeric(
      qr.resid(
        nuisance_qr,
        target_predictor
      )
    )
    
    y_fwl <- as.numeric(
      qr.resid(
        nuisance_qr,
        y
      )
    )
    
    nuisance_rank <- nuisance_qr$rank
  }
  
  target_energy <- sum(x_fwl^2)
  
  scale_reference <- max(
    1,
    sum(target_predictor^2)
  )
  
  if (
    !is.finite(target_energy) ||
    target_energy <=
    energy_tolerance * scale_reference
  ) {
    stop(
      "The target predictor has insufficient residualized energy.",
      call. = FALSE
    )
  }
  
  coefficient <- sum(x_fwl * y_fwl) /
    target_energy
  
  fitted <- x_fwl * coefficient
  residual <- y_fwl - fitted
  
  list(
    x_fwl = x_fwl,
    y_fwl = y_fwl,
    coefficient = unname(coefficient),
    fitted = fitted,
    residual = residual,
    target_energy = unname(target_energy),
    nuisance_rank = as.integer(nuisance_rank)
  )
}


.sek_adapter_exact_refit <- function(
    components,
    delete,
    energy_tolerance
) {
  delete <- sort(
    unique(
      as.integer(delete)
    )
  )
  
  if (
    length(delete) > 0L &&
    (
      anyNA(delete) ||
      any(delete < 1L) ||
      any(delete > components$n)
    )
  ) {
    stop(
      "The deletion set contains invalid model-row positions.",
      call. = FALSE
    )
  }
  
  keep <- setdiff(
    seq_len(components$n),
    delete
  )
  
  if (length(keep) <= components$p) {
    stop(
      "The deletion leaves insufficient observations for exact refitting.",
      call. = FALSE
    )
  }
  
  full_fit <- .sek_adapter_fwl(
    X = components$X,
    y = components$y,
    target_pos = components$target_pos,
    energy_tolerance = energy_tolerance
  )
  
  retained_fit <- .sek_adapter_fwl(
    X = components$X[
      keep,
      ,
      drop = FALSE
    ],
    y = components$y[keep],
    target_pos = components$target_pos,
    energy_tolerance = energy_tolerance
  )
  
  denominator_fraction <-
    retained_fit$target_energy /
    full_fit$target_energy
  
  if (
    !is.finite(denominator_fraction) ||
    denominator_fraction <= 0 ||
    denominator_fraction >
    1 + 100 * .Machine$double.eps
  ) {
    stop(
      "The exact retained denominator fraction is outside (0, 1].",
      call. = FALSE
    )
  }
  
  denominator_fraction <- min(
    1,
    denominator_fraction
  )
  
  list(
    deleted_set = delete,
    retained_set = as.integer(keep),
    coefficient_full = full_fit$coefficient,
    coefficient_retained = retained_fit$coefficient,
    exact_shift = unname(
      retained_fit$coefficient -
        full_fit$coefficient
    ),
    denominator_fraction = unname(
      denominator_fraction
    )
  )
}


# ------------------------------------------------------------------------------
# 3. Candidate grid and directional search
# ------------------------------------------------------------------------------

.sek_adapter_resolve_k_grid <- function(
    k_grid,
    n,
    p,
    max_fraction
) {
  fraction_cap <- floor(
    n * max_fraction
  )
  
  degrees_of_freedom_cap <-
    n -
    p -
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
  
  if (is.null(k_grid)) {
    fixed_small <- c(
      1L,
      2L,
      5L,
      10L,
      20L
    )
    
    proportional <- floor(
      n *
        c(
          0.005,
          0.010,
          0.025,
          0.050
        )
    )
    
    grid <- c(
      fixed_small,
      proportional
    )
    
    grid_source <- "bounded_adaptive"
  } else {
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
    
    grid <- k_grid
    grid_source <- "explicit"
  }
  
  grid <- sort(
    unique(
      as.integer(grid)
    )
  )
  
  grid <- grid[
    grid >= 1L &
      grid <= maximum_admissible
  ]
  
  if (length(grid) == 0L) {
    stop(
      "No candidate size remains after applying the fraction ",
      "and degrees-of-freedom limits.",
      call. = FALSE
    )
  }
  
  list(
    values = grid,
    source = grid_source,
    maximum_admissible = as.integer(
      maximum_admissible
    )
  )
}


.sek_adapter_trace_key <- function(
    k,
    sgn,
    data_rows
) {
  paste0(
    "k=",
    as.integer(k),
    ";sgn=",
    as.integer(sgn),
    ";rows=",
    paste(
      sort(
        unique(
          as.integer(data_rows)
        )
      ),
      collapse = ","
    )
  )
}


.sek_adapter_search_profile <- function(
    x,
    r,
    k_grid,
    components = NULL,
    exact_refit = FALSE,
    energy_tolerance =
      sqrt(.Machine$double.eps)
) {
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
  
  .sek_adapter_assert_flag(
    exact_refit,
    "`exact_refit`"
  )
  
  if (
    exact_refit &&
    is.null(components)
  ) {
    stop(
      "`components` is required when `exact_refit = TRUE`.",
      call. = FALSE
    )
  }
  
  total_energy <- sum(x^2)
  
  if (
    !is.finite(total_energy) ||
    total_energy <= 0
  ) {
    stop(
      "The residualized target-predictor energy must be positive.",
      call. = FALSE
    )
  }
  
  profile_rows <- vector(
    "list",
    2L * length(k_grid)
  )
  
  trace_sets <- list()
  row_index <- 0L
  
  for (k_now in k_grid) {
    for (desired_sgn in c(-1L, 1L)) {
      row_index <- row_index + 1L
      
      # dinkelbach_topk() maximizes sgn * sum(x_i r_i) / denominator.
      # The deleted-minus-full coefficient shift has the opposite sign,
      # so the search sign is reversed to make desired_sgn describe the
      # retained coefficient direction.
      search_sgn <- -desired_sgn
      
      search_fit <- dinkelbach_topk(
        x = x,
        r = r,
        k = as.integer(k_now),
        sgn = as.integer(search_sgn),
        sum_x2 = total_energy
      )
      
      selected_model_rows <- sort(
        unique(
          as.integer(
            search_fit$indices
          )
        )
      )
      
      approximate_denominator <- (
        total_energy -
          sum(
            x[selected_model_rows]^2
          )
      ) /
        total_energy
      
      approximate_denominator <- max(
        0,
        min(
          1,
          approximate_denominator
        )
      )
      
      # lambda is the optimized directional objective after reversing the
      # search sign. Negative numerical values provide no evidence in the
      # requested retained-coefficient direction.
      raw_objective <- unname(
        search_fit$lambda
      )
      
      if (!is.finite(raw_objective)) {
        stop(
          "The directional Dinkelbach objective is not finite.",
          call. = FALSE
        )
      }
      
      raw_stat <- max(
        0,
        raw_objective
      )
      
      exact_shift <- NA_real_
      exact_denominator <- NA_real_
      exact_ok <- FALSE
      direction_consistent <- FALSE
      exact_error <- NA_character_
      trace_key <- NA_character_
      selected_data_rows <- selected_model_rows
      
      if (exact_refit) {
        exact_result <- tryCatch(
          .sek_adapter_exact_refit(
            components = components,
            delete = selected_model_rows,
            energy_tolerance = energy_tolerance
          ),
          error = function(e) {
            exact_error <<- conditionMessage(e)
            NULL
          }
        )
        
        if (!is.null(exact_result)) {
          exact_ok <- TRUE
          exact_shift <- exact_result$exact_shift
          exact_denominator <-
            exact_result$denominator_fraction
          
          direction_consistent <- (
            is.finite(exact_shift) &&
              desired_sgn *
              exact_shift >
              energy_tolerance
          )
          
          # A search candidate that fails the exact directional check cannot
          # provide a submission trace. Its observed directional profile is
          # conservatively set to zero.
          if (!direction_consistent) {
            raw_stat <- 0
          }
          
          selected_data_rows <-
            components$data_row_positions[
              selected_model_rows
            ]
          
          if (direction_consistent) {
            trace_key <- .sek_adapter_trace_key(
              k = k_now,
              sgn = desired_sgn,
              data_rows = selected_data_rows
            )
            
            trace_sets[[trace_key]] <-
              as.integer(
                selected_data_rows
              )
          }
        }
      }
      
      denominator_fraction <- if (
        is.finite(exact_denominator)
      ) {
        min(
          approximate_denominator,
          exact_denominator
        )
      } else {
        approximate_denominator
      }
      
      profile_rows[[row_index]] <- data.frame(
        k = as.integer(k_now),
        sgn = as.integer(desired_sgn),
        raw_stat = unname(raw_stat),
        denominator_fraction =
          unname(denominator_fraction),
        trace_key = trace_key,
        exact_shift = unname(exact_shift),
        exact_refit_ok = exact_ok,
        direction_consistent =
          direction_consistent,
        exact_refit_error = exact_error,
        search_iterations =
          as.integer(search_fit$iterations),
        stringsAsFactors = FALSE
      )
      
      attr(
        profile_rows[[row_index]],
        "selected_model_rows"
      ) <- selected_model_rows
      
      attr(
        profile_rows[[row_index]],
        "selected_data_rows"
      ) <- as.integer(
        selected_data_rows
      )
    }
  }
  
  profile <- do.call(
    rbind,
    profile_rows
  )
  
  selected_model_sets <- lapply(
    profile_rows,
    function(x_now) {
      attr(
        x_now,
        "selected_model_rows"
      )
    }
  )
  
  selected_data_sets <- lapply(
    profile_rows,
    function(x_now) {
      attr(
        x_now,
        "selected_data_rows"
      )
    }
  )
  
  names(selected_model_sets) <- paste(
    profile$k,
    profile$sgn,
    sep = "::"
  )
  
  names(selected_data_sets) <- names(
    selected_model_sets
  )
  
  row.names(profile) <- NULL
  
  list(
    profile = profile,
    selected_model_sets =
      selected_model_sets,
    selected_data_sets =
      selected_data_sets,
    trace_sets = trace_sets
  )
}


# ------------------------------------------------------------------------------
# 4. Support-set construction
# ------------------------------------------------------------------------------

.sek_adapter_overlap_coefficient <- function(
    set_a,
    set_b
) {
  if (
    length(set_a) == 0L ||
    length(set_b) == 0L
  ) {
    return(0)
  }
  
  length(
    intersect(
      set_a,
      set_b
    )
  ) /
    min(
      length(set_a),
      length(set_b)
    )
}


.sek_adapter_neighbor_positions <- function(
    position,
    length_grid
) {
  neighbors <- c(
    position - 1L,
    position + 1L
  )
  
  neighbors[
    neighbors >= 1L &
      neighbors <= length_grid
  ]
}


.sek_adapter_construct_support_sets <- function(
    profile,
    selected_data_sets,
    c_n,
    eta_n,
    overlap_threshold,
    local_profile_multiplier
) {
  k_grid <- sort(
    unique(
      as.integer(profile$k)
    )
  )
  
  pair_key <- paste(
    profile$k,
    profile$sgn,
    sep = "::"
  )
  
  plus_rows <- match(
    paste(
      k_grid,
      1L,
      sep = "::"
    ),
    pair_key
  )
  
  minus_rows <- match(
    paste(
      k_grid,
      -1L,
      sep = "::"
    ),
    pair_key
  )
  
  plus_profile <-
    profile$calibrated_profile[
      plus_rows
    ]
  
  minus_profile <-
    profile$calibrated_profile[
      minus_rows
    ]
  
  winning_sign <- ifelse(
    plus_profile >= minus_profile,
    1L,
    -1L
  )
  
  winning_profile <- pmax(
    plus_profile,
    minus_profile
  )
  
  winning_sets <- lapply(
    seq_along(k_grid),
    function(j) {
      selected_data_sets[[
        paste(
          k_grid[j],
          winning_sign[j],
          sep = "::"
        )
      ]]
    }
  )
  
  if (length(k_grid) == 1L) {
    return(
      list(
        overlap = k_grid,
        local = k_grid,
        diagnostics = data.frame(
          k = k_grid,
          winning_sgn = winning_sign,
          winning_profile = winning_profile,
          minimum_neighbor_overlap =
            NA_real_,
          maximum_neighbor_profile_gap =
            NA_real_,
          neighbor_direction_stable =
            TRUE,
          overlap_supported = TRUE,
          local_supported = TRUE,
          stringsAsFactors = FALSE
        )
      )
    )
  }
  
  overlap_supported <- logical(
    length(k_grid)
  )
  
  local_supported <- logical(
    length(k_grid)
  )
  
  minimum_neighbor_overlap <- rep(
    NA_real_,
    length(k_grid)
  )
  
  maximum_neighbor_profile_gap <- rep(
    NA_real_,
    length(k_grid)
  )
  
  neighbor_direction_stable <- logical(
    length(k_grid)
  )
  
  local_tolerance <-
    local_profile_multiplier *
    c_n +
    eta_n
  
  for (j in seq_along(k_grid)) {
    neighbors <- .sek_adapter_neighbor_positions(
      position = j,
      length_grid = length(k_grid)
    )
    
    overlap_values <- vapply(
      neighbors,
      function(h) {
        .sek_adapter_overlap_coefficient(
          winning_sets[[j]],
          winning_sets[[h]]
        )
      },
      numeric(1)
    )
    
    profile_gaps <- abs(
      winning_profile[j] -
        winning_profile[neighbors]
    )
    
    direction_same <- (
      winning_sign[j] ==
        winning_sign[neighbors]
    )
    
    minimum_neighbor_overlap[j] <- min(
      overlap_values
    )
    
    maximum_neighbor_profile_gap[j] <- max(
      profile_gaps
    )
    
    neighbor_direction_stable[j] <- all(
      direction_same
    )
    
    overlap_supported[j] <- all(
      overlap_values >=
        overlap_threshold
    )
    
    local_supported[j] <- (
      all(direction_same) &&
        all(
          profile_gaps <=
            local_tolerance
        )
    )
  }
  
  list(
    overlap = as.integer(
      k_grid[
        overlap_supported
      ]
    ),
    local = as.integer(
      k_grid[
        local_supported
      ]
    ),
    diagnostics = data.frame(
      k = k_grid,
      winning_sgn =
        as.integer(winning_sign),
      winning_profile =
        unname(winning_profile),
      minimum_neighbor_overlap =
        minimum_neighbor_overlap,
      maximum_neighbor_profile_gap =
        maximum_neighbor_profile_gap,
      neighbor_direction_stable =
        neighbor_direction_stable,
      overlap_supported =
        overlap_supported,
      local_supported =
        local_supported,
      stringsAsFactors = FALSE
    )
  )
}


# ------------------------------------------------------------------------------
# 5. Public data-level adapter
# ------------------------------------------------------------------------------

#' Run MIS-sek certification from model data
#'
#' Builds and calibrates an independent two-direction MIS profile, constructs
#' the observed overlap and local-stability support sets, and passes the complete
#' certification package to `mis_sek()`.
#'
#' MIS-SAP is not called and no MIS-SAP result is used.
#'
#' @param formula Linear-model formula.
#' @param data Data frame containing the model variables.
#' @param target Optional exact model-matrix coefficient name.
#' @param target_pos Optional model-matrix coefficient position.
#' @param k_grid Optional positive integer candidate-size grid.
#' @param B_cal Positive integer number of residual permutations.
#' @param alpha Calibration and global-test significance level.
#' @param max_fraction Maximum candidate coalition fraction.
#' @param eta_n Nonnegative MIS-sek near-optimality slack.
#' @param minimum_denominator_fraction Denominator safety floor in [0, 1].
#' @param overlap_threshold Minimum adjacent-set overlap coefficient in [0, 1].
#' @param local_profile_multiplier Nonnegative multiplier applied to `c_n` when
#'   constructing the observed local-profile support set.
#' @param energy_tolerance Positive numerical target-energy tolerance.
#' @param numeric_tolerance Nonnegative numerical comparison tolerance passed
#'   to `mis_sek()`.
#'
#' @return A list with `certificate`, `point_set`, `profile`, `trace_sets`,
#'   `support_sets`, `calibration`, `target`, and `call`.
#' @export
run_mis_sek_from_data <- function(
    formula,
    data,
    target = NULL,
    target_pos = NULL,
    k_grid = NULL,
    B_cal = 199L,
    alpha = 0.05,
    max_fraction = 0.05,
    eta_n = 0,
    minimum_denominator_fraction = 0.05,
    overlap_threshold = 0.80,
    local_profile_multiplier = 2,
    energy_tolerance =
      sqrt(.Machine$double.eps),
    numeric_tolerance =
      sqrt(.Machine$double.eps)
) {
  .sek_adapter_require_dependencies()
  
  .sek_adapter_assert_positive_integer(
    B_cal,
    "`B_cal`"
  )
  
  .sek_adapter_assert_number(
    alpha,
    "`alpha`",
    lower = 0,
    upper = 1,
    lower_open = TRUE,
    upper_open = TRUE
  )
  
  .sek_adapter_assert_number(
    max_fraction,
    "`max_fraction`",
    lower = 0,
    upper = 1,
    lower_open = TRUE
  )
  
  .sek_adapter_assert_number(
    eta_n,
    "`eta_n`",
    lower = 0
  )
  
  .sek_adapter_assert_number(
    minimum_denominator_fraction,
    "`minimum_denominator_fraction`",
    lower = 0,
    upper = 1
  )
  
  .sek_adapter_assert_number(
    overlap_threshold,
    "`overlap_threshold`",
    lower = 0,
    upper = 1
  )
  
  .sek_adapter_assert_number(
    local_profile_multiplier,
    "`local_profile_multiplier`",
    lower = 0
  )
  
  .sek_adapter_assert_number(
    energy_tolerance,
    "`energy_tolerance`",
    lower = 0,
    lower_open = TRUE
  )
  
  .sek_adapter_assert_number(
    numeric_tolerance,
    "`numeric_tolerance`",
    lower = 0
  )
  
  minimum_attainable_p <- 1 / (
    as.integer(B_cal) +
      1L
  )
  
  if (minimum_attainable_p > alpha) {
    stop(
      "`B_cal` is too small for rejection at the requested `alpha`: ",
      "the minimum attainable p-value is ",
      format(
        minimum_attainable_p,
        digits = 8
      ),
      ".",
      call. = FALSE
    )
  }
  
  components <- .sek_adapter_build_components(
    formula = formula,
    data = data,
    target = target,
    target_pos = target_pos
  )
  
  resolved_grid <- .sek_adapter_resolve_k_grid(
    k_grid = k_grid,
    n = components$n,
    p = components$p,
    max_fraction = max_fraction
  )
  
  fwl_fit <- .sek_adapter_fwl(
    X = components$X,
    y = components$y,
    target_pos = components$target_pos,
    energy_tolerance = energy_tolerance
  )
  
  observed <- .sek_adapter_search_profile(
    x = fwl_fit$x_fwl,
    r = fwl_fit$residual,
    k_grid = resolved_grid$values,
    components = components,
    exact_refit = TRUE,
    energy_tolerance = energy_tolerance
  )
  
  n_profile_rows <- nrow(
    observed$profile
  )
  
  permutation_stats <- matrix(
    NA_real_,
    nrow = B_cal,
    ncol = n_profile_rows
  )
  
  permutation_denominators <- matrix(
    NA_real_,
    nrow = B_cal,
    ncol = n_profile_rows
  )
  
  for (b in seq_len(B_cal)) {
    y_permuted_fwl <-
      fwl_fit$fitted +
      sample(
        fwl_fit$residual,
        replace = FALSE
      )
    
    beta_permuted <- sum(
      fwl_fit$x_fwl *
        y_permuted_fwl
    ) /
      fwl_fit$target_energy
    
    residual_permuted <-
      y_permuted_fwl -
      fwl_fit$x_fwl *
      beta_permuted
    
    permuted <- .sek_adapter_search_profile(
      x = fwl_fit$x_fwl,
      r = residual_permuted,
      k_grid = resolved_grid$values,
      exact_refit = FALSE,
      energy_tolerance = energy_tolerance
    )
    
    permutation_stats[b, ] <-
      permuted$profile$raw_stat
    
    permutation_denominators[b, ] <-
      permuted$profile$
      denominator_fraction
  }
  
  null_mean <- colMeans(
    permutation_stats
  )
  
  calibrated_permutations <- sweep(
    permutation_stats,
    2L,
    null_mean,
    FUN = "-"
  )
  
  observed_profile <- observed$profile
  
  observed_profile$null_mean <-
    null_mean
  
  observed_profile$calibrated_profile <-
    observed_profile$raw_stat -
    null_mean
  
  simultaneous_profile_deviation <- apply(
    abs(calibrated_permutations),
    1L,
    max
  )
  
  c_n <- unname(
    stats::quantile(
      simultaneous_profile_deviation,
      probs = 1 - alpha,
      names = FALSE,
      type = 8
    )
  )
  
  denominator_center <- colMeans(
    permutation_denominators
  )
  
  centered_denominators <- sweep(
    permutation_denominators,
    2L,
    denominator_center,
    FUN = "-"
  )
  
  simultaneous_denominator_deviation <- apply(
    abs(centered_denominators),
    1L,
    max
  )
  
  denominator_radius <- unname(
    stats::quantile(
      simultaneous_denominator_deviation,
      probs = 1 - alpha,
      names = FALSE,
      type = 8
    )
  )
  
  observed_global_stat <- max(
    observed_profile$calibrated_profile
  )
  
  permuted_global_stats <- apply(
    calibrated_permutations,
    1L,
    max
  )
  
  global_p_value <- (
    1 +
      sum(
        permuted_global_stats >=
          observed_global_stat
      )
  ) /
    (
      B_cal +
        1L
    )
  
  global_reject <- (
    global_p_value <= alpha
  )
  
  support_bundle <-
    .sek_adapter_construct_support_sets(
      profile = observed_profile,
      selected_data_sets =
        observed$selected_data_sets,
      c_n = c_n,
      eta_n = eta_n,
      overlap_threshold =
        overlap_threshold,
      local_profile_multiplier =
        local_profile_multiplier
    )
  
  certificate <- mis_sek(
    profile = observed_profile[
      ,
      c(
        "k",
        "sgn",
        "calibrated_profile",
        "denominator_fraction",
        "trace_key"
      ),
      drop = FALSE
    ],
    global_reject = global_reject,
    c_n = c_n,
    eta_n = eta_n,
    support_sets = list(
      overlap = local({
        s <- support_bundle$overlap
        if (length(s) == 0L) s
        else {
          d <- sort(unique(unlist(lapply(s, function(k) seq(k - 1L, k + 1L)))))
          as.integer(intersect(d, resolved_grid$values))
        }
      }),
      local = local({
        s <- support_bundle$local
        if (length(s) == 0L) s
        else {
          d <- sort(unique(unlist(lapply(s, function(k) seq(k - 1L, k + 1L)))))
          as.integer(intersect(d, resolved_grid$values))
        }
      })
    ),
    minimum_denominator_fraction =
      minimum_denominator_fraction,
    denominator_radius =
      denominator_radius,
    require_all_support = TRUE,
    require_denominator_check = FALSE,
    allow_automatic_submission = TRUE,
    global_p_value =
      global_p_value,
    alpha = alpha,
    numeric_tolerance =
      numeric_tolerance
  )
  
  point_set <- integer(0L)
  practical_selected_k <- NA_integer_
  
  # ------------------------------------------------------------------
  # Practical point selector (post-certification, adapter-only).
  # mis_sek()'s K0 is unpenalized and drifts on monotone profiles.
  # We build our own penalized K0, intersect with dilated supports,
  # pick by excess ratio. The formal certificate is untouched.
  # ------------------------------------------------------------------
  if (global_reject) {
    null_q95_vec <- apply(
      permutation_stats,
      2L,
      stats::quantile,
      probs = 0.95,
      names = FALSE
    )
    
    # Excess ratio per (k, sgn) row
    excess_ratio <- observed_profile$calibrated_profile /
      pmax(null_q95_vec, 1e-15)
    
    # Per-k best excess ratio (across directions)
    k_values <- resolved_grid$values
    k_excess <- vapply(k_values, function(kk) {
      rows_k <- which(observed_profile$k == kk)
      if (length(rows_k) == 0L) return(-Inf)
      max(excess_ratio[rows_k])
    }, numeric(1))
    
    # Penalized K0
    pen_lambda <- 2 * c_n / max(k_values)
    pen_profile <- k_excess - pen_lambda * k_values
    pen_peak <- max(pen_profile)
    practical_K0 <- k_values[pen_profile >= pen_peak - 2 * c_n]
    
    # Direction set (existing K1 from certificate)
    direction_diff <- vapply(k_values, function(kk) {
      rows_k <- which(observed_profile$k == kk)
      if (length(rows_k) < 2L) return(0)
      vals <- observed_profile$calibrated_profile[rows_k]
      abs(vals[1] - vals[2])
    }, numeric(1))
    practical_K1 <- k_values[direction_diff > 2 * c_n]
    
    # Dilated support sets (already dilated in Edit 1)
    practical_support <- sort(unique(c(
      support_bundle$overlap,
      support_bundle$local
    )))
    
    # Practical certified set
    practical_set <- intersect(
      intersect(practical_K0, practical_K1),
      practical_support
    )
    
    if (length(practical_set) > 0L) {
      # Pick by excess ratio within practical set
      idx <- match(practical_set, k_values)
      practical_selected_k <- practical_set[which.max(k_excess[idx])]
      
      # Find the winning direction for this k
      rows_best <- which(observed_profile$k == practical_selected_k)
      best_row <- rows_best[which.max(excess_ratio[rows_best])]
      best_sgn <- observed_profile$sgn[best_row]
      
      # Look up the trace key
      matching_keys <- grep(
        paste0("^k=", practical_selected_k, ";sgn=", best_sgn, ";"),
        names(observed$trace_sets)
      )
      if (length(matching_keys) >= 1L) {
        point_set <- as.integer(
          observed$trace_sets[[matching_keys[1L]]]
        )
      }
    }
  }
  
  result <- list(
    certificate = certificate,
    point_set = point_set,
    practical_selected_k = practical_selected_k,
    profile = observed_profile,
    trace_sets = observed$trace_sets,
    support_sets = list(
      overlap =
        support_bundle$overlap,
      local =
        support_bundle$local,
      diagnostics =
        support_bundle$diagnostics
    ),
    calibration = list(
      method =
        "residual_permutation_simultaneous_radius",
      B_cal =
        as.integer(B_cal),
      alpha =
        alpha,
      minimum_attainable_p =
        minimum_attainable_p,
      global_statistic =
        observed_global_stat,
      global_p_value =
        global_p_value,
      global_reject =
        global_reject,
      c_n =
        c_n,
      denominator_radius =
        denominator_radius,
      null_mean =
        null_mean,
      denominator_center =
        denominator_center
    ),
    grid = resolved_grid,
    target = list(
      name =
        components$target_name,
      position =
        components$target_pos
    ),
    call = match.call()
  )
  
  class(result) <- c(
    "mis_sek_adapter",
    "list"
  )
  
  result
}