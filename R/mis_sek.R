# ============================================================
# MIS-sek — Set-valued stable effective-k selector
#
# Four support sets:
#
#   K0 = calibrated-profile eligible set
#   K1 = direction-consistency support set
#   K2 = overlap-stability support set
#   K3 = local-smoothness support set
#
# Jointly supported set:
#
#   K_MIS-sek = K0 intersect K1 intersect K2 intersect K3
#
# The support sets do not vote and are not averaged.
#
# Possible states remain:
#
#   - "no_excessive_influence"
#   - "stable_effective_k"
#   - "excessive_influence_but_k_unstable"
#
# A stable result may be:
#
#   - point identified: one jointly supported k;
#   - set identified: several jointly supported k values.
#
# Automatic top-k submission is allowed only for a singleton.
# ============================================================


#' Select a stable effective influence scale
#'
#' @param profile Data frame containing the calibrated MIS profile.
#'   It must contain k, calibrated-profile, and sign columns.
#' @param global_reject Logical scalar indicating whether the global
#'   excessive-influence test rejected.
#' @param c_n Nonnegative stochastic calibration tolerance.
#' @param eta_n Nonnegative deterministic slack tolerance.
#' @param support_sets Named list containing the three supporting
#'   sets: direction, overlap, and local.
#' @param minimum_denominator_fraction Minimum acceptable denominator
#'   fraction.
#' @param require_all_support Whether all three support sets must be
#'   supplied and nonempty.
#' @param require_denominator_check Whether denominator safety must
#'   be verified.
#' @param k_col Name of the k column.
#' @param profile_col Name of the calibrated-profile column.
#' @param sgn_col Name of the direction column.
#' @param denominator_col Optional denominator-fraction column.
#' @param trace_key_col Optional trace-key column.
#' @param global_p_value Optional global p-value for reporting.
#' @param alpha Optional global significance level.
#' @param tie_tolerance Numerical tolerance for ties.
#'
#' @return An object of class "mis_sek".
#'
#' @export
mis_sek <- function(
    profile,
    global_reject,
    c_n,
    eta_n,
    support_sets = list(
      direction = NULL,
      overlap = NULL,
      local = NULL
    ),
    minimum_denominator_fraction = 0.05,
    require_all_support = TRUE,
    require_denominator_check = TRUE,
    k_col = "k",
    profile_col = "calibrated_profile",
    sgn_col = "sgn",
    denominator_col = "denominator_fraction",
    trace_key_col = "trace_key",
    global_p_value = NA_real_,
    alpha = NA_real_,
    tie_tolerance = sqrt(.Machine$double.eps)
) {
  # ----------------------------------------------------------
  # Validate basic inputs
  # ----------------------------------------------------------
  if (
    !is.data.frame(profile) ||
    nrow(profile) == 0L
  ) {
    stop(
      "profile must be a nonempty data frame."
    )
  }
  
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
      )
    )
  }
  
  if (
    length(global_reject) != 1L ||
    is.na(global_reject)
  ) {
    stop(
      "global_reject must be one non-missing logical value."
    )
  }
  
  global_reject <- as.logical(
    global_reject
  )
  
  if (
    length(c_n) != 1L ||
    !is.finite(c_n) ||
    c_n < 0
  ) {
    stop(
      "c_n must be one finite nonnegative number."
    )
  }
  
  if (
    length(eta_n) != 1L ||
    !is.finite(eta_n) ||
    eta_n < 0
  ) {
    stop(
      "eta_n must be one finite nonnegative number."
    )
  }
  
  if (
    length(minimum_denominator_fraction) != 1L ||
    !is.finite(minimum_denominator_fraction) ||
    minimum_denominator_fraction < 0 ||
    minimum_denominator_fraction > 1
  ) {
    stop(
      "minimum_denominator_fraction must lie between 0 and 1."
    )
  }
  
  if (!is.list(support_sets)) {
    stop(
      "support_sets must be a named list."
    )
  }
  
  
  # ----------------------------------------------------------
  # Validate profile values
  # ----------------------------------------------------------
  profile_k <- as.numeric(
    profile[[k_col]]
  )
  
  profile_value <- as.numeric(
    profile[[profile_col]]
  )
  
  profile_sgn <- as.numeric(
    profile[[sgn_col]]
  )
  
  if (
    any(!is.finite(profile_k)) ||
    any(profile_k < 1) ||
    any(
      abs(
        profile_k -
        round(profile_k)
      ) > tie_tolerance
    )
  ) {
    stop(
      "All k values must be positive finite integers."
    )
  }
  
  if (any(!is.finite(profile_value))) {
    stop(
      "All calibrated profile values must be finite."
    )
  }
  
  if (
    any(!is.finite(profile_sgn)) ||
    any(
      !profile_sgn %in%
      c(-1, 1)
    )
  ) {
    stop(
      "The sign column must contain only -1 and 1."
    )
  }
  
  profile_k <- as.integer(
    round(profile_k)
  )
  
  profile_sgn <- as.integer(
    profile_sgn
  )
  
  k_grid <- sort(
    unique(profile_k)
  )
  
  
  # ----------------------------------------------------------
  # Optional denominator information
  # ----------------------------------------------------------
  denominator_available <- (
    !is.null(denominator_col) &&
      denominator_col %in%
      names(profile)
  )
  
  if (denominator_available) {
    profile_denominator <- as.numeric(
      profile[[denominator_col]]
    )
    
    if (
      any(!is.finite(profile_denominator)) ||
      any(profile_denominator < 0)
    ) {
      stop(
        "Denominator fractions must be finite and nonnegative."
      )
    }
  } else {
    profile_denominator <- rep(
      NA_real_,
      nrow(profile)
    )
  }
  
  
  # ----------------------------------------------------------
  # Collapse both directions to one profile value per k
  # ----------------------------------------------------------
  collapsed_rows <- lapply(
    k_grid,
    function(k_now) {
      positions <- which(
        profile_k == k_now
      )
      
      values_now <- profile_value[
        positions
      ]
      
      maximum_now <- max(
        values_now
      )
      
      tie_band <- (
        tie_tolerance *
          max(
            1,
            abs(maximum_now)
          )
      )
      
      tied_positions <- positions[
        abs(
          values_now -
            maximum_now
        ) <= tie_band
      ]
      
      tied_signs <- sort(
        unique(
          profile_sgn[
            tied_positions
          ]
        )
      )
      
      direction_unique <- (
        length(tied_signs) == 1L
      )
      
      winning_sign <- if (
        direction_unique
      ) {
        tied_signs[1L]
      } else {
        NA_integer_
      }
      
      # Conservative denominator value when tied.
      denominator_now <- if (
        denominator_available
      ) {
        min(
          profile_denominator[
            tied_positions
          ]
        )
      } else {
        NA_real_
      }
      
      trace_key_now <- NA_character_
      
      if (
        trace_key_col %in%
        names(profile)
      ) {
        trace_key_now <- as.character(
          profile[
            tied_positions[1L],
            trace_key_col
          ]
        )
      }
      
      data.frame(
        k =
          k_now,
        
        calibrated_profile =
          maximum_now,
        
        sgn =
          winning_sign,
        
        direction_unique =
          direction_unique,
        
        denominator_fraction =
          denominator_now,
        
        trace_key =
          trace_key_now,
        
        number_of_tied_rows =
          length(tied_positions),
        
        number_of_tied_directions =
          length(tied_signs),
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  collapsed_profile <- do.call(
    rbind,
    collapsed_rows
  )
  
  rownames(collapsed_profile) <- NULL
  
  
  # ----------------------------------------------------------
  # K0: calibrated-profile eligible support set
  # ----------------------------------------------------------
  selection_tolerance <- (
    2 * c_n +
      eta_n
  )
  
  maximum_profile <- max(
    collapsed_profile$
      calibrated_profile
  )
  
  profile_cutoff <- (
    maximum_profile -
      selection_tolerance
  )
  
  cutoff_tolerance <- (
    tie_tolerance *
      max(
        1,
        abs(maximum_profile)
      )
  )
  
  primary_support_set <- sort(
    collapsed_profile$k[
      collapsed_profile$
        calibrated_profile >=
        profile_cutoff -
        cutoff_tolerance
    ]
  )
  
  primary_reference_k <- min(
    primary_support_set
  )
  
  
  # ----------------------------------------------------------
  # Normalize K1, K2, and K3
  # ----------------------------------------------------------
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
        "An unnamed support_sets list must have length 3."
      )
    }
    
    names(support_sets) <-
      expected_support_names
  }
  
  unknown_support_names <- setdiff(
    names(support_sets),
    expected_support_names
  )
  
  if (
    length(unknown_support_names) > 0L
  ) {
    stop(
      "Unknown support-set names: ",
      paste(
        unknown_support_names,
        collapse = ", "
      )
    )
  }
  
  normalized_support_sets <- setNames(
    vector(
      "list",
      length(expected_support_names)
    ),
    expected_support_names
  )
  
  normalize_k_set <- function(
    x,
    name
  ) {
    if (
      is.null(x) ||
      length(x) == 0L
    ) {
      return(integer(0L))
    }
    
    x <- as.numeric(x)
    
    if (
      any(!is.finite(x)) ||
      any(x < 1) ||
      any(
        abs(
          x -
          round(x)
        ) > tie_tolerance
      )
    ) {
      stop(
        name,
        " support set must contain positive finite integers."
      )
    }
    
    x <- sort(
      unique(
        as.integer(
          round(x)
        )
      )
    )
    
    outside_grid <- setdiff(
      x,
      k_grid
    )
    
    if (length(outside_grid) > 0L) {
      stop(
        name,
        " support set contains values outside the k grid: ",
        paste(
          outside_grid,
          collapse = ", "
        )
      )
    }
    
    x
  }
  
  for (
    support_name in
    expected_support_names
  ) {
    normalized_support_sets[[
      support_name
    ]] <- normalize_k_set(
      support_sets[[
        support_name
      ]],
      support_name
    )
  }
  
  support_available <- vapply(
    normalized_support_sets,
    length,
    integer(1L)
  ) > 0L
  
  support_complete <- all(
    support_available
  )
  
  
  # ----------------------------------------------------------
  # Joint support-set intersection
  # ----------------------------------------------------------
  if (
    require_all_support &&
    !support_complete
  ) {
    identified_k_set <- integer(0L)
  } else {
    active_support_sets <-
      normalized_support_sets[
        support_available
      ]
    
    sets_to_intersect <- c(
      list(
        primary =
          primary_support_set
      ),
      active_support_sets
    )
    
    identified_k_set <- sort(
      Reduce(
        intersect,
        sets_to_intersect
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Diagnostic union
  #
  # Used only for reporting when the intersection is empty.
  # It is not an identified set or confidence set.
  # ----------------------------------------------------------
  diagnostic_candidate_k_set <- sort(
    unique(
      unlist(
        c(
          list(
            primary_support_set
          ),
          normalized_support_sets
        ),
        use.names = FALSE
      )
    )
  )
  
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
  
  identified_k_range <- make_range(
    identified_k_set
  )
  
  diagnostic_candidate_k_range <-
    make_range(
      diagnostic_candidate_k_set
    )
  
  
  # ----------------------------------------------------------
  # Direction consistency across the intersection
  # ----------------------------------------------------------
  identified_positions <- match(
    identified_k_set,
    collapsed_profile$k
  )
  
  direction_stable <- FALSE
  identified_sgn <- NA_integer_
  
  if (length(identified_k_set) > 0L) {
    identified_direction_unique <- all(
      collapsed_profile$
        direction_unique[
          identified_positions
        ]
    )
    
    identified_signs <- unique(
      collapsed_profile$sgn[
        identified_positions
      ]
    )
    
    identified_signs <- identified_signs[
      is.finite(identified_signs)
    ]
    
    direction_stable <- (
      identified_direction_unique &&
        length(identified_signs) == 1L
    )
    
    if (direction_stable) {
      identified_sgn <- as.integer(
        identified_signs[1L]
      )
    }
  }
  
  
  # ----------------------------------------------------------
  # Denominator safety across the intersection
  # ----------------------------------------------------------
  denominator_safe <- TRUE
  minimum_observed_denominator <- NA_real_
  
  if (require_denominator_check) {
    if (
      !denominator_available ||
      length(identified_k_set) == 0L
    ) {
      denominator_safe <- FALSE
    } else {
      denominator_values <- collapsed_profile$
        denominator_fraction[
          identified_positions
        ]
      
      minimum_observed_denominator <- min(
        denominator_values
      )
      
      denominator_safe <- all(
        denominator_values >=
          minimum_denominator_fraction
      )
    }
  } else if (
    denominator_available &&
    length(identified_k_set) > 0L
  ) {
    minimum_observed_denominator <- min(
      collapsed_profile$
        denominator_fraction[
          identified_positions
        ]
    )
  }
  
  
  # ----------------------------------------------------------
  # Stability reasons
  # ----------------------------------------------------------
  instability_reasons <- character(0L)
  
  if (
    require_all_support &&
    !support_complete
  ) {
    instability_reasons <- c(
      instability_reasons,
      "missing_support_set"
    )
  }
  
  if (
    global_reject &&
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
      "direction_not_stable_over_intersection"
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
  
  
  # ----------------------------------------------------------
  # Final three-state decision
  # ----------------------------------------------------------
  if (!global_reject) {
    state <- "no_excessive_influence"
    
    selection_type <- "none"
    
    selected_k <- NA_integer_
    selected_k_set <- integer(0L)
    selected_sgn <- NA_integer_
    selected_trace_key <- NA_character_
    
    automatic_top_k_submission <- FALSE
  } else if (
    length(identified_k_set) > 0L &&
    direction_stable &&
    denominator_safe
  ) {
    state <- "stable_effective_k"
    
    selected_k_set <- identified_k_set
    
    if (length(identified_k_set) == 1L) {
      selection_type <- "point"
      
      selected_k <- identified_k_set[1L]
      selected_sgn <- identified_sgn
      
      selected_position <- match(
        selected_k,
        collapsed_profile$k
      )
      
      selected_trace_key <- collapsed_profile$
        trace_key[
          selected_position
        ]
      
      automatic_top_k_submission <- TRUE
    } else {
      selection_type <- "set"
      
      selected_k <- NA_integer_
      selected_sgn <- NA_integer_
      selected_trace_key <- NA_character_
      
      automatic_top_k_submission <- FALSE
    }
  } else {
    state <-
      "excessive_influence_but_k_unstable"
    
    selection_type <- "none"
    
    selected_k <- NA_integer_
    selected_k_set <- integer(0L)
    selected_sgn <- NA_integer_
    selected_trace_key <- NA_character_
    
    automatic_top_k_submission <- FALSE
  }
  
  
  # ----------------------------------------------------------
  # Membership table
  # ----------------------------------------------------------
  support_membership <- data.frame(
    k =
      k_grid,
    
    primary = (
      k_grid %in%
        primary_support_set
    ),
    
    direction = if (
      support_available["direction"]
    ) {
      k_grid %in%
        normalized_support_sets$
        direction
    } else {
      NA
    },
    
    overlap = if (
      support_available["overlap"]
    ) {
      k_grid %in%
        normalized_support_sets$
        overlap
    } else {
      NA
    },
    
    local = if (
      support_available["local"]
    ) {
      k_grid %in%
        normalized_support_sets$
        local
    } else {
      NA
    },
    
    jointly_supported = (
      k_grid %in%
        identified_k_set
    ),
    
    stringsAsFactors = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Result
  # ----------------------------------------------------------
  result <- list(
    method =
      "MIS-sek",
    
    state =
      state,
    
    selection_type =
      selection_type,
    
    selected_k =
      selected_k,
    
    selected_k_set =
      selected_k_set,
    
    selected_k_range =
      make_range(
        selected_k_set
      ),
    
    selected_sgn =
      selected_sgn,
    
    selected_trace_key =
      selected_trace_key,
    
    automatic_top_k_submission =
      automatic_top_k_submission,
    
    primary = list(
      reference_k =
        primary_reference_k,
      
      support_set =
        primary_support_set,
      
      support_range =
        make_range(
          primary_support_set
        ),
      
      maximum_profile =
        maximum_profile,
      
      profile_cutoff =
        profile_cutoff,
      
      c_n =
        c_n,
      
      eta_n =
        eta_n,
      
      selection_tolerance =
        selection_tolerance
    ),
    
    support_sets =
      normalized_support_sets,
    
    identified = list(
      k_set =
        identified_k_set,
      
      k_range =
        identified_k_range,
      
      sgn =
        identified_sgn,
      
      direction_stable =
        direction_stable
    ),
    
    diagnostic = list(
      candidate_k_set =
        diagnostic_candidate_k_set,
      
      candidate_k_range =
        diagnostic_candidate_k_range,
      
      warning = paste(
        "The diagnostic candidate range is a bounding range",
        "of the reported support sets, not a confidence interval."
      )
    ),
    
    global = list(
      reject =
        global_reject,
      
      p_value =
        global_p_value,
      
      alpha =
        alpha
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
      
      support_available =
        support_available,
      
      support_complete =
        support_complete,
      
      direction_stable =
        direction_stable,
      
      denominator_available =
        denominator_available,
      
      denominator_safe =
        denominator_safe,
      
      minimum_observed_denominator =
        minimum_observed_denominator,
      
      required_minimum_denominator =
        minimum_denominator_fraction,
      
      instability_reasons =
        instability_reasons
    ),
    
    support_membership =
      support_membership,
    
    collapsed_profile =
      collapsed_profile,
    
    call =
      match.call()
  )
  
  class(result) <- c(
    "mis_sek",
    "list"
  )
  
  result
}


# ------------------------------------------------------------
# Print method
# ------------------------------------------------------------
#' @export
print.mis_sek <- function(
    x,
    ...
) {
  format_k_set <- function(k_set) {
    if (length(k_set) == 0L) {
      return("<empty>")
    }
    
    paste(
      k_set,
      collapse = ", "
    )
  }
  
  format_k_range <- function(k_range) {
    if (anyNA(k_range)) {
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
  
  cat("MIS-sek set-valued effective-k selector\n")
  cat("---------------------------------------\n")
  cat("State:", x$state, "\n")
  
  if (is.finite(x$global$p_value)) {
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
    "Primary support set:",
    format_k_set(
      x$primary$support_set
    ),
    "\n"
  )
  
  cat(
    "Direction support set:",
    format_k_set(
      x$support_sets$direction
    ),
    "\n"
  )
  
  cat(
    "Overlap support set:",
    format_k_set(
      x$support_sets$overlap
    ),
    "\n"
  )
  
  cat(
    "Local support set:",
    format_k_set(
      x$support_sets$local
    ),
    "\n"
  )
  
  cat(
    "Joint intersection:",
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
      "No effective k is selected.\n"
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
      "Selected direction:",
      x$selected_sgn,
      "\n"
    )
    
    cat(
      "Automatic top-k submission: yes\n"
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
      "Stable possible-k set:",
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
      "Automatic top-k submission: no\n"
    )
  } else {
    cat(
      "No jointly supported effective k.\n"
    )
    
    cat(
      "Diagnostic candidate values:",
      format_k_set(
        x$diagnostic$
          candidate_k_set
      ),
      "\n"
    )
    
    cat(
      "Diagnostic bounding range:",
      format_k_range(
        x$diagnostic$
          candidate_k_range
      ),
      "\n"
    )
    
    cat(
      "Instability reasons:",
      paste(
        x$stability$
          instability_reasons,
        collapse = ", "
      ),
      "\n"
    )
    
    cat(
      "Automatic top-k submission: no\n"
    )
  }
  
  invisible(x)
}