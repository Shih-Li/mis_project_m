# ==============================================================================
# MIS-sek: theorem-aligned conservative effective-k certification
#
# Replace the existing mis_sek() and print.mis_sek() definitions with this file.
#
# Core contract:
#   K_sek_hat = K0_hat ∩ K1_hat ∩ K2_hat ∩ K3_hat
#
#   K0_hat: calibrated near-optimal profile set, constructed only on the
#           denominator-admissible domain.
#   K1_hat: internally certified direction set:
#           |M_hat(k,+1) - M_hat(k,-1)| > 2*c_n.
#   K2_hat: observed overlap support set.
#   K3_hat: observed local-stability support set.
#
# A computed empty support set remains empty and annihilates the intersection.
# NULL means "not computed" and produces diagnostics_incomplete.
# ==============================================================================


#' MIS-sek theorem-aligned effective-k certification
#'
#' @param profile Nonempty data frame containing exactly one row for every
#'   `(k, sgn)` pair, with both `sgn = -1` and `sgn = 1` for every k.
#' @param global_reject One non-missing logical value.
#' @param c_n Finite nonnegative simultaneous profile error radius.
#' @param eta_n Finite nonnegative near-optimality slack.
#' @param support_sets Named list with `overlap` and `local`. A legacy
#'   `direction` component may be supplied; when non-NULL, it must exactly
#'   equal the internally certified direction set.
#' @param minimum_denominator_fraction Denominator safety threshold in [0, 1].
#' @param denominator_radius Finite nonnegative simultaneous denominator
#'   error radius. A row is empirically admissible only when
#'   `denominator >= minimum_denominator_fraction + denominator_radius`.
#' @param require_all_support Must be TRUE in the theorem-aligned interface.
#' @param require_denominator_check Whether denominator admissibility is
#'   enforced before constructing K0 and K1.
#' @param allow_automatic_submission Whether a fully certified singleton with
#'   a non-missing trace key may be automatically submitted.
#' @param k_col,profile_col,sgn_col,denominator_col,trace_key_col Column names.
#' @param global_p_value Optional scalar p-value in [0, 1].
#' @param alpha Optional scalar level in [0, 1].
#' @param numeric_tolerance Floating-point comparison tolerance. It is not a
#'   statistical uncertainty radius.
#' @param tie_tolerance Deprecated legacy argument. When supplied, it is used
#'   only as `numeric_tolerance`.
#'
#' @return Object of class `mis_sek`.
#' @export
mis_sek <- function(
    profile,
    global_reject,
    c_n,
    eta_n = 0,
    support_sets = list(
      direction = NULL,
      overlap = NULL,
      local = NULL
    ),
    minimum_denominator_fraction = 0.05,
    denominator_radius = 0,
    require_all_support = TRUE,
    require_denominator_check = TRUE,
    allow_automatic_submission = FALSE,
    k_col = "k",
    profile_col = "calibrated_profile",
    sgn_col = "sgn",
    denominator_col = "denominator_fraction",
    trace_key_col = "trace_key",
    global_p_value = NA_real_,
    alpha = NA_real_,
    numeric_tolerance = sqrt(.Machine$double.eps),
    tie_tolerance = NULL
) {
  # ---------------------------------------------------------------------------
  # 0. Strict scalar helpers
  # ---------------------------------------------------------------------------
  assert_flag <- function(x, name) {
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
  
  assert_number <- function(
    x,
    name,
    lower = -Inf,
    upper = Inf,
    allow_na = FALSE
  ) {
    if (
      allow_na &&
      length(x) == 1L &&
      is.numeric(x) &&
      is.na(x)
    ) {
      return(invisible(TRUE))
    }
    
    if (
      !is.numeric(x) ||
      length(x) != 1L ||
      !is.finite(x) ||
      x < lower ||
      x > upper
    ) {
      stop(
        name,
        " must be one finite numeric value in [",
        lower,
        ", ",
        upper,
        "].",
        call. = FALSE
      )
    }
    
    invisible(TRUE)
  }
  
  assert_numeric_column <- function(x, name) {
    if (
      is.factor(x) ||
      is.character(x) ||
      is.logical(x) ||
      !is.numeric(x)
    ) {
      stop(
        name,
        " must be stored as a numeric/integer column; ",
        "silent factor or character coercion is not allowed.",
        call. = FALSE
      )
    }
    
    invisible(TRUE)
  }
  
  make_range <- function(x) {
    if (length(x) == 0L) {
      return(
        c(
          lower = NA_integer_,
          upper = NA_integer_
        )
      )
    }
    
    c(
      lower = min(x),
      upper = max(x)
    )
  }
  
  normalize_integer_set <- function(
    x,
    name,
    k_grid,
    numeric_tolerance
  ) {
    if (is.null(x)) {
      return(NULL)
    }
    
    if (length(x) == 0L) {
      return(integer(0L))
    }
    
    if (
      is.factor(x) ||
      is.character(x) ||
      is.logical(x) ||
      !is.numeric(x)
    ) {
      stop(
        name,
        " must be a numeric/integer vector or integer(0); ",
        "factor and character coercion is not allowed.",
        call. = FALSE
      )
    }
    
    if (
      any(!is.finite(x)) ||
      any(x < 1) ||
      any(
        abs(
          x -
          round(x)
        ) >
        numeric_tolerance
      )
    ) {
      stop(
        name,
        " must contain positive finite integers.",
        call. = FALSE
      )
    }
    
    output <- sort(
      unique(
        as.integer(
          round(x)
        )
      )
    )
    
    outside <- setdiff(
      output,
      k_grid
    )
    
    if (length(outside) > 0L) {
      stop(
        name,
        " contains values outside the profile k-grid: ",
        paste(
          outside,
          collapse = ", "
        ),
        call. = FALSE
      )
    }
    
    output
  }
  
  
  # ---------------------------------------------------------------------------
  # 1. Strict scalar and interface validation
  # ---------------------------------------------------------------------------
  if (
    !is.data.frame(profile) ||
    nrow(profile) == 0L
  ) {
    stop(
      "profile must be a nonempty data frame.",
      call. = FALSE
    )
  }
  
  assert_flag(
    global_reject,
    "global_reject"
  )
  
  assert_flag(
    require_all_support,
    "require_all_support"
  )
  
  assert_flag(
    require_denominator_check,
    "require_denominator_check"
  )
  
  assert_flag(
    allow_automatic_submission,
    "allow_automatic_submission"
  )
  
  if (!require_all_support) {
    stop(
      "The formal MIS-sek interface requires require_all_support = TRUE. ",
      "Partial diagnostic intersections do not have the formal guarantee.",
      call. = FALSE
    )
  }
  
  assert_number(
    c_n,
    "c_n",
    lower = 0
  )
  
  assert_number(
    eta_n,
    "eta_n",
    lower = 0
  )
  
  assert_number(
    minimum_denominator_fraction,
    "minimum_denominator_fraction",
    lower = 0,
    upper = 1
  )
  
  assert_number(
    denominator_radius,
    "denominator_radius",
    lower = 0
  )
  
  assert_number(
    numeric_tolerance,
    "numeric_tolerance",
    lower = 0
  )
  
  assert_number(
    global_p_value,
    "global_p_value",
    lower = 0,
    upper = 1,
    allow_na = TRUE
  )
  
  assert_number(
    alpha,
    "alpha",
    lower = 0,
    upper = 1,
    allow_na = TRUE
  )
  
  if (!is.null(tie_tolerance)) {
    assert_number(
      tie_tolerance,
      "tie_tolerance",
      lower = 0
    )
    
    warning(
      "tie_tolerance is deprecated. In MIS-sek it is used only ",
      "as numeric_tolerance and never as a statistical radius.",
      call. = FALSE
    )
    
    numeric_tolerance <- tie_tolerance
  }
  
  if (
    is.finite(global_p_value) &&
    is.finite(alpha)
  ) {
    implied_reject <- (
      global_p_value <= alpha
    )
    
    if (!identical(global_reject, implied_reject)) {
      stop(
        "global_reject is inconsistent with global_p_value <= alpha.",
        call. = FALSE
      )
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # 2. Validate required columns and strict profile schema
  # ---------------------------------------------------------------------------
  required_columns <- c(
    k_col,
    profile_col,
    sgn_col
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(profile)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "profile is missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  assert_numeric_column(
    profile[[k_col]],
    k_col
  )
  
  assert_numeric_column(
    profile[[profile_col]],
    profile_col
  )
  
  assert_numeric_column(
    profile[[sgn_col]],
    sgn_col
  )
  
  profile_k_raw <- profile[[k_col]]
  profile_value <- as.numeric(
    profile[[profile_col]]
  )
  profile_sgn_raw <- profile[[sgn_col]]
  
  if (
    any(!is.finite(profile_k_raw)) ||
    any(profile_k_raw < 1) ||
    any(
      abs(
        profile_k_raw -
        round(profile_k_raw)
      ) >
      numeric_tolerance
    )
  ) {
    stop(
      "All k values must be positive finite integers.",
      call. = FALSE
    )
  }
  
  if (any(!is.finite(profile_value))) {
    stop(
      "All calibrated profile values must be finite.",
      call. = FALSE
    )
  }
  
  if (
    any(!is.finite(profile_sgn_raw)) ||
    any(
      !profile_sgn_raw %in%
      c(-1, 1)
    )
  ) {
    stop(
      "The sign column must contain only -1 and 1.",
      call. = FALSE
    )
  }
  
  profile_k <- as.integer(
    round(profile_k_raw)
  )
  
  profile_sgn <- as.integer(
    profile_sgn_raw
  )
  
  profile_pair_key <- paste(
    profile_k,
    profile_sgn,
    sep = "::"
  )
  
  duplicated_pairs <- unique(
    profile_pair_key[
      duplicated(profile_pair_key)
    ]
  )
  
  if (length(duplicated_pairs) > 0L) {
    stop(
      "profile must contain exactly one row per (k, sgn) pair. ",
      "Duplicated pairs: ",
      paste(
        duplicated_pairs,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  k_grid <- sort(
    unique(profile_k)
  )
  
  expected_pair_keys <- as.vector(
    outer(
      k_grid,
      c(-1L, 1L),
      function(k_now, sgn_now) {
        paste(
          k_now,
          sgn_now,
          sep = "::"
        )
      }
    )
  )
  
  missing_pair_keys <- setdiff(
    expected_pair_keys,
    profile_pair_key
  )
  
  if (length(missing_pair_keys) > 0L) {
    stop(
      "Every k must have exactly one sgn = -1 row and one sgn = 1 row. ",
      "Missing pairs: ",
      paste(
        missing_pair_keys,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  plus_rows <- match(
    paste(
      k_grid,
      1L,
      sep = "::"
    ),
    profile_pair_key
  )
  
  minus_rows <- match(
    paste(
      k_grid,
      -1L,
      sep = "::"
    ),
    profile_pair_key
  )
  
  
  # ---------------------------------------------------------------------------
  # 3. Denominator validation and admissible domain
  #
  # This occurs BEFORE K0 construction.
  # ---------------------------------------------------------------------------
  denominator_available <- (
    !is.null(denominator_col) &&
      length(denominator_col) == 1L &&
      denominator_col %in%
      names(profile)
  )
  
  if (
    require_denominator_check &&
    !denominator_available
  ) {
    stop(
      "require_denominator_check = TRUE, but denominator_col is unavailable.",
      call. = FALSE
    )
  }
  
  if (denominator_available) {
    assert_numeric_column(
      profile[[denominator_col]],
      denominator_col
    )
    
    profile_denominator <- as.numeric(
      profile[[denominator_col]]
    )
    
    if (
      any(!is.finite(profile_denominator)) ||
      any(profile_denominator < 0) ||
      any(profile_denominator > 1)
    ) {
      stop(
        "Denominator fractions must be finite values in [0, 1].",
        call. = FALSE
      )
    }
  } else {
    profile_denominator <- rep(
      NA_real_,
      nrow(profile)
    )
  }
  
  empirical_denominator_cutoff <- (
    minimum_denominator_fraction +
      denominator_radius
  )
  
  if (
    require_denominator_check &&
    empirical_denominator_cutoff > 1
  ) {
    stop(
      "minimum_denominator_fraction + denominator_radius exceeds 1; ",
      "no row can be certified denominator-safe.",
      call. = FALSE
    )
  }
  
  row_admissible <- if (
    require_denominator_check
  ) {
    profile_denominator >=
      empirical_denominator_cutoff
  } else {
    rep(
      TRUE,
      nrow(profile)
    )
  }
  
  # Conservative theorem-aligned choice: both directions must be denominator
  # admissible before k can enter K0 or K1.
  k_admissible <- (
    row_admissible[plus_rows] &
      row_admissible[minus_rows]
  )
  
  
  # ---------------------------------------------------------------------------
  # 4. Exact two-direction profile and internally certified K1
  # ---------------------------------------------------------------------------
  profile_plus <- profile_value[
    plus_rows
  ]
  
  profile_minus <- profile_value[
    minus_rows
  ]
  
  profile_max <- pmax(
    profile_plus,
    profile_minus
  )
  
  signed_direction_difference <- (
    profile_plus -
      profile_minus
  )
  
  direction_gap <- abs(
    signed_direction_difference
  )
  
  direction_certification_threshold <- (
    2 *
      c_n
  )
  
  direction_certified <- (
    k_admissible &
      direction_gap >
      direction_certification_threshold +
      numeric_tolerance
  )
  
  winning_sign <- ifelse(
    signed_direction_difference > 0,
    1L,
    -1L
  )
  
  certified_sign <- ifelse(
    direction_certified,
    winning_sign,
    NA_integer_
  )
  
  denominator_plus <- profile_denominator[
    plus_rows
  ]
  
  denominator_minus <- profile_denominator[
    minus_rows
  ]
  
  minimum_denominator_by_k <- if (
    denominator_available
  ) {
    pmin(
      denominator_plus,
      denominator_minus
    )
  } else {
    rep(
      NA_real_,
      length(k_grid)
    )
  }
  
  trace_available <- (
    !is.null(trace_key_col) &&
      length(trace_key_col) == 1L &&
      trace_key_col %in%
      names(profile)
  )
  
  trace_values <- if (
    trace_available
  ) {
    as.character(
      profile[[trace_key_col]]
    )
  } else {
    rep(
      NA_character_,
      nrow(profile)
    )
  }
  
  trace_plus <- trace_values[
    plus_rows
  ]
  
  trace_minus <- trace_values[
    minus_rows
  ]
  
  selected_direction_trace <- ifelse(
    winning_sign == 1L,
    trace_plus,
    trace_minus
  )
  
  collapsed_profile <- data.frame(
    k = k_grid,
    profile_plus = profile_plus,
    profile_minus = profile_minus,
    calibrated_profile = profile_max,
    direction_gap = direction_gap,
    direction_certified = direction_certified,
    sgn = certified_sign,
    denominator_plus = denominator_plus,
    denominator_minus = denominator_minus,
    denominator_fraction = minimum_denominator_by_k,
    denominator_admissible = k_admissible,
    trace_key = selected_direction_trace,
    stringsAsFactors = FALSE
  )
  
  
  # ---------------------------------------------------------------------------
  # 5. K0: calibrated near-optimal set on the admissible domain
  # ---------------------------------------------------------------------------
  selection_tolerance <- (
    2 *
      c_n +
      eta_n
  )
  
  if (any(k_admissible)) {
    maximum_profile <- max(
      profile_max[
        k_admissible
      ]
    )
    
    profile_cutoff <- (
      maximum_profile -
        selection_tolerance
    )
    
    cutoff_numeric_tolerance <- (
      numeric_tolerance *
        max(
          1,
          abs(maximum_profile)
        )
    )
    
    primary_support_set <- sort(
      k_grid[
        k_admissible &
          profile_max >=
          profile_cutoff -
          cutoff_numeric_tolerance
      ]
    )
  } else {
    maximum_profile <- NA_real_
    profile_cutoff <- NA_real_
    cutoff_numeric_tolerance <- numeric_tolerance
    primary_support_set <- integer(0L)
  }
  
  primary_reference_k <- if (
    length(primary_support_set) > 0L
  ) {
    min(primary_support_set)
  } else {
    NA_integer_
  }
  
  direction_support_set <- sort(
    k_grid[
      direction_certified
    ]
  )
  
  
  # ---------------------------------------------------------------------------
  # 6. Normalize K2 and K3 while preserving NULL versus integer(0)
  # ---------------------------------------------------------------------------
  if (!is.list(support_sets)) {
    stop(
      "support_sets must be a named list.",
      call. = FALSE
    )
  }
  
  expected_support_names <- c(
    "direction",
    "overlap",
    "local"
  )
  
  if (is.null(names(support_sets))) {
    if (
      length(support_sets) !=
      length(expected_support_names)
    ) {
      stop(
        "An unnamed support_sets list must have length 3.",
        call. = FALSE
      )
    }
    
    names(support_sets) <- expected_support_names
  }
  
  if (
    anyNA(names(support_sets)) ||
    any(names(support_sets) == "")
  ) {
    stop(
      "Every support_sets component must have a nonempty name.",
      call. = FALSE
    )
  }
  
  duplicated_support_names <- unique(
    names(support_sets)[
      duplicated(
        names(support_sets)
      )
    ]
  )
  
  if (
    length(duplicated_support_names) >
    0L
  ) {
    stop(
      "Duplicated support-set names are not allowed: ",
      paste(
        duplicated_support_names,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  unknown_support_names <- setdiff(
    names(support_sets),
    expected_support_names
  )
  
  if (
    length(unknown_support_names) >
    0L
  ) {
    stop(
      "Unknown support-set names: ",
      paste(
        unknown_support_names,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  for (
    support_name in
    setdiff(
      expected_support_names,
      names(support_sets)
    )
  ) {
    support_sets[[
      support_name
    ]] <- NULL
  }
  
  legacy_direction_set <- normalize_integer_set(
    support_sets$direction,
    "direction support set",
    k_grid,
    numeric_tolerance
  )
  
  if (
    !is.null(legacy_direction_set) &&
    !identical(
      legacy_direction_set,
      direction_support_set
    )
  ) {
    stop(
      "The supplied direction support set disagrees with the internally ",
      "certified set defined by |M_hat(k,+1)-M_hat(k,-1)| > 2*c_n. ",
      "Remove support_sets$direction or update its upstream construction.",
      call. = FALSE
    )
  }
  
  overlap_support_set <- normalize_integer_set(
    support_sets$overlap,
    "overlap support set",
    k_grid,
    numeric_tolerance
  )
  
  local_support_set <- normalize_integer_set(
    support_sets$local,
    "local support set",
    k_grid,
    numeric_tolerance
  )
  
  support_available <- c(
    direction = TRUE,
    overlap = !is.null(
      overlap_support_set
    ),
    local = !is.null(
      local_support_set
    )
  )
  
  support_complete <- all(
    support_available
  )
  
  normalized_support_sets <- list(
    direction = direction_support_set,
    overlap = overlap_support_set,
    local = local_support_set
  )
  
  
  # ---------------------------------------------------------------------------
  # 7. Exact hard intersection
  #
  # Computed empty sets are included. Missing sets are not silently dropped.
  # ---------------------------------------------------------------------------
  if (support_complete) {
    identified_k_set <- sort(
      Reduce(
        intersect,
        list(
          primary_support_set,
          direction_support_set,
          overlap_support_set,
          local_support_set
        )
      )
    )
  } else {
    identified_k_set <- integer(0L)
  }
  
  identified_k_range <- make_range(
    identified_k_set
  )
  
  identified_positions <- match(
    identified_k_set,
    collapsed_profile$k
  )
  
  identified_signs <- if (
    length(identified_k_set) > 0L
  ) {
    unique(
      collapsed_profile$sgn[
        identified_positions
      ]
    )
  } else {
    integer(0L)
  }
  
  identified_signs <- identified_signs[
    !is.na(identified_signs)
  ]
  
  direction_stable <- (
    length(identified_k_set) > 0L &&
      length(identified_signs) == 1L
  )
  
  identified_sgn <- if (
    direction_stable
  ) {
    as.integer(
      identified_signs[1L]
    )
  } else {
    NA_integer_
  }
  
  minimum_observed_denominator <- if (
    denominator_available &&
    length(identified_k_set) > 0L
  ) {
    min(
      collapsed_profile$denominator_fraction[
        identified_positions
      ]
    )
  } else {
    NA_real_
  }
  
  denominator_safe <- if (
    require_denominator_check
  ) {
    length(identified_k_set) > 0L &&
      all(
        collapsed_profile$denominator_admissible[
          identified_positions
        ]
      )
  } else {
    TRUE
  }
  
  
  # ---------------------------------------------------------------------------
  # 8. State classification and reasons
  # ---------------------------------------------------------------------------
  instability_reasons <- character(0L)
  
  if (!support_complete) {
    instability_reasons <- c(
      instability_reasons,
      "missing_support_set"
    )
  }
  
  if (
    require_denominator_check &&
    !any(k_admissible)
  ) {
    instability_reasons <- c(
      instability_reasons,
      "no_denominator_admissible_k"
    )
  }
  
  if (
    global_reject &&
    support_complete &&
    length(identified_k_set) == 0L
  ) {
    instability_reasons <- c(
      instability_reasons,
      "empty_support_intersection"
    )
  }
  
  if (
    global_reject &&
    length(identified_k_set) > 0L &&
    !direction_stable
  ) {
    instability_reasons <- c(
      instability_reasons,
      "direction_conflict_over_intersection"
    )
  }
  
  if (
    global_reject &&
    length(identified_k_set) > 0L &&
    !denominator_safe
  ) {
    instability_reasons <- c(
      instability_reasons,
      "denominator_not_safe"
    )
  }
  
  instability_reasons <- unique(
    instability_reasons
  )
  
  state <- "no_stable_effective_k"
  selection_type <- "none"
  selected_k <- NA_integer_
  selected_k_set <- integer(0L)
  selected_sgn <- NA_integer_
  selected_trace_key <- NA_character_
  automatic_top_k_submission <- FALSE
  
  if (!global_reject) {
    state <- "no_excessive_influence"
  } else if (!support_complete) {
    state <- "diagnostics_incomplete"
  } else if (
    length(identified_k_set) > 0L &&
    direction_stable &&
    denominator_safe
  ) {
    state <- "stable_effective_k"
    selected_k_set <- identified_k_set
    selected_sgn <- identified_sgn
    
    if (length(identified_k_set) == 1L) {
      selection_type <- "point"
      selected_k <- identified_k_set[1L]
      
      selected_position <- match(
        selected_k,
        collapsed_profile$k
      )
      
      selected_trace_key <- collapsed_profile$trace_key[
        selected_position
      ]
      
      trace_is_usable <- (
        length(selected_trace_key) == 1L &&
          !is.na(selected_trace_key) &&
          nzchar(selected_trace_key)
      )
      
      automatic_top_k_submission <- (
        allow_automatic_submission &&
          trace_is_usable
      )
    } else {
      selection_type <- "set"
      # Keep the common certified direction for set-valued output.
      selected_k <- NA_integer_
      selected_trace_key <- NA_character_
      automatic_top_k_submission <- FALSE
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # 9. Reporting tables and result object
  # ---------------------------------------------------------------------------
  support_membership <- data.frame(
    k = k_grid,
    primary = (
      k_grid %in%
        primary_support_set
    ),
    direction = (
      k_grid %in%
        direction_support_set
    ),
    overlap = if (
      is.null(overlap_support_set)
    ) {
      NA
    } else {
      k_grid %in%
        overlap_support_set
    },
    local = if (
      is.null(local_support_set)
    ) {
      NA
    } else {
      k_grid %in%
        local_support_set
    },
    jointly_supported = if (
      support_complete
    ) {
      k_grid %in%
        identified_k_set
    } else {
      NA
    },
    stringsAsFactors = FALSE
  )
  
  diagnostic_candidate_k_set <- sort(
    unique(
      c(
        primary_support_set,
        direction_support_set,
        if (
          is.null(overlap_support_set)
        ) {
          integer(0L)
        } else {
          overlap_support_set
        },
        if (
          is.null(local_support_set)
        ) {
          integer(0L)
        } else {
          local_support_set
        }
      )
    )
  )
  
  formal_guarantee_eligible <- (
    support_complete &&
      require_denominator_check &&
      denominator_radius >= 0 &&
      c_n >= 0
  )
  
  result <- list(
    method = "MIS-sek",
    state = state,
    selection_type = selection_type,
    selected_k = selected_k,
    selected_k_set = selected_k_set,
    selected_k_range = make_range(
      selected_k_set
    ),
    selected_sgn = selected_sgn,
    selected_trace_key = selected_trace_key,
    automatic_top_k_submission = automatic_top_k_submission,
    
    primary = list(
      reference_k = primary_reference_k,
      support_set = primary_support_set,
      support_range = make_range(
        primary_support_set
      ),
      maximum_profile = maximum_profile,
      profile_cutoff = profile_cutoff,
      c_n = c_n,
      eta_n = eta_n,
      selection_tolerance = selection_tolerance,
      admissible_k_set = k_grid[
        k_admissible
      ]
    ),
    
    direction = list(
      support_set = direction_support_set,
      threshold = direction_certification_threshold,
      rule = "abs(profile_plus - profile_minus) > 2*c_n"
    ),
    
    support_sets = normalized_support_sets,
    
    identified = list(
      k_set = identified_k_set,
      k_range = identified_k_range,
      sgn = identified_sgn,
      direction_stable = direction_stable
    ),
    
    diagnostic = list(
      candidate_k_set = diagnostic_candidate_k_set,
      candidate_k_range = make_range(
        diagnostic_candidate_k_set
      ),
      warning = paste(
        "The diagnostic candidate range is a bounding range",
        "of reported support sets, not a confidence interval."
      )
    ),
    
    global = list(
      reject = global_reject,
      p_value = global_p_value,
      alpha = alpha
    ),
    
    stability = list(
      stable = (
        state ==
          "stable_effective_k"
      ),
      point_identified = (
        state ==
          "stable_effective_k" &&
          selection_type ==
          "point"
      ),
      set_identified = (
        state ==
          "stable_effective_k" &&
          selection_type ==
          "set"
      ),
      support_available = support_available,
      support_complete = support_complete,
      direction_stable = direction_stable,
      denominator_available = denominator_available,
      denominator_safe = denominator_safe,
      denominator_radius = denominator_radius,
      empirical_denominator_cutoff = empirical_denominator_cutoff,
      minimum_observed_denominator = minimum_observed_denominator,
      required_minimum_denominator = minimum_denominator_fraction,
      formal_guarantee_eligible = formal_guarantee_eligible,
      instability_reasons = instability_reasons
    ),
    
    support_membership = support_membership,
    collapsed_profile = collapsed_profile,
    call = match.call()
  )
  
  class(result) <- c(
    "mis_sek",
    "list"
  )
  
  result
}


#' @export
print.mis_sek <- function(
    x,
    ...
) {
  format_k_set <- function(k_set) {
    if (is.null(k_set)) {
      return("<not computed>")
    }
    
    if (length(k_set) == 0L) {
      return("<empty>")
    }
    
    paste(
      k_set,
      collapse = ", "
    )
  }
  
  format_k_range <- function(k_range) {
    if (
      is.null(k_range) ||
      anyNA(k_range)
    ) {
      return("<undefined>")
    }
    
    paste0(
      "[",
      k_range["lower"],
      ", ",
      k_range["upper"],
      "]"
    )
  }
  
  cat("MIS-sek theorem-aligned effective-k certification\n")
  cat("----------------------------------------------------\n")
  cat("State:", x$state, "\n")
  
  if (
    is.numeric(x$global$p_value) &&
    length(x$global$p_value) == 1L &&
    is.finite(x$global$p_value)
  ) {
    cat(
      "Global p-value:",
      format(
        x$global$p_value,
        digits = 4
      ),
      "\n"
    )
  }
  
  cat(
    "K0 primary support:",
    format_k_set(
      x$primary$support_set
    ),
    "\n"
  )
  
  cat(
    "K1 certified direction support:",
    format_k_set(
      x$direction$support_set
    ),
    "\n"
  )
  
  cat(
    "K2 overlap support:",
    format_k_set(
      x$support_sets$overlap
    ),
    "\n"
  )
  
  cat(
    "K3 local support:",
    format_k_set(
      x$support_sets$local
    ),
    "\n"
  )
  
  cat(
    "Hard intersection:",
    format_k_set(
      x$identified$k_set
    ),
    "\n"
  )
  
  if (
    identical(
      x$state,
      "no_excessive_influence"
    )
  ) {
    cat(
      "The global test did not reject; no effective k is selected.\n"
    )
  } else if (
    identical(
      x$state,
      "diagnostics_incomplete"
    )
  ) {
    cat(
      "The overlap or local diagnostic was not computed; ",
      "no formal MIS-sek decision is available.\n",
      sep = ""
    )
  } else if (
    identical(
      x$state,
      "stable_effective_k"
    ) &&
    identical(
      x$selection_type,
      "point"
    )
  ) {
    cat(
      "Selected effective k:",
      x$selected_k,
      "\n"
    )
    
    cat(
      "Certified direction:",
      x$selected_sgn,
      "\n"
    )
    
    cat(
      "Automatic top-k submission:",
      if (
        isTRUE(
          x$automatic_top_k_submission
        )
      ) {
        "yes"
      } else {
        "no"
      },
      "\n"
    )
  } else if (
    identical(
      x$state,
      "stable_effective_k"
    ) &&
    identical(
      x$selection_type,
      "set"
    )
  ) {
    cat(
      "Stable effective-k set:",
      format_k_set(
        x$selected_k_set
      ),
      "\n"
    )
    
    cat(
      "Bounding range:",
      format_k_range(
        x$selected_k_range
      ),
      "\n"
    )
    
    cat(
      "Common certified direction:",
      x$selected_sgn,
      "\n"
    )
    
    cat(
      "Automatic top-k submission: no\n"
    )
  } else {
    cat(
      "No stable effective k was certified.\n"
    )
    
    if (
      length(
        x$stability$instability_reasons
      ) > 0L
    ) {
      cat(
        "Reasons:",
        paste(
          x$stability$instability_reasons,
          collapse = ", "
        ),
        "\n"
      )
    }
    
    cat(
      "Diagnostic candidate values:",
      format_k_set(
        x$diagnostic$candidate_k_set
      ),
      "\n"
    )
    
    cat(
      "Automatic top-k submission: no\n"
    )
  }
  
  invisible(x)
}