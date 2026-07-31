# ==============================================================================
# File: /R/test-mis_sek.R
# Purpose: Regression and contract tests for the MIS-sek effective-k
#          certification procedure. Verifies strict profile validation,
#          denominator-first admissibility, statistically certified direction
#          support, exact hard-intersection semantics, explicit abstention
#          states, stable provenance, and monotonicity under support-set
#          tightening.
#
# Dependencies: Requires /R/mis_sek.R and the testthat package.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Test fixtures
# ------------------------------------------------------------------------------

make_valid_profile <- function(
    k_grid = 1:3,
    plus = c(1.0, 0.9, 0.8),
    minus = c(0.2, 0.2, 0.2),
    denominator = 0.5
) {
  stopifnot(
    length(k_grid) == length(plus),
    length(k_grid) == length(minus)
  )
  
  data.frame(
    k = rep(
      k_grid,
      each = 2L
    ),
    sgn = rep(
      c(-1L, 1L),
      times = length(k_grid)
    ),
    calibrated_profile = as.vector(
      rbind(
        minus,
        plus
      )
    ),
    denominator_fraction = denominator,
    trace_key = paste0(
      "row-",
      seq_len(
        2L *
          length(k_grid)
      )
    ),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 1. Support-set completeness and empty-set semantics
# ------------------------------------------------------------------------------

testthat::test_that(
  "computed empty support annihilates the intersection",
  {
    profile <- make_valid_profile()
    
    result <- mis_sek(
      profile = profile,
      global_reject = TRUE,
      c_n = 0.01,
      eta_n = 0,
      support_sets = list(
        overlap = integer(0L),
        local = 1:3
      ),
      denominator_radius = 0.01
    )
    
    testthat::expect_identical(
      result$identified$k_set,
      integer(0L)
    )
    
    testthat::expect_identical(
      result$state,
      "no_stable_effective_k"
    )
  }
)


testthat::test_that(
  "missing support produces diagnostics_incomplete",
  {
    profile <- make_valid_profile()
    
    result <- mis_sek(
      profile = profile,
      global_reject = TRUE,
      c_n = 0.01,
      support_sets = list(
        overlap = NULL,
        local = 1:3
      ),
      denominator_radius = 0.01
    )
    
    testthat::expect_identical(
      result$state,
      "diagnostics_incomplete"
    )
  }
)

# ------------------------------------------------------------------------------
# 2. Profile schema and input validation
# ------------------------------------------------------------------------------

testthat::test_that(
  "both directions are required for every k",
  {
    profile <- make_valid_profile()
    profile <- profile[-1L, ]
    
    testthat::expect_error(
      mis_sek(
        profile = profile,
        global_reject = TRUE,
        c_n = 0.01,
        support_sets = list(
          overlap = 1:3,
          local = 1:3
        )
      ),
      "Every k must have exactly one"
    )
  }
)


testthat::test_that(
  "duplicate k-sign rows are rejected",
  {
    profile <- make_valid_profile()
    profile <- rbind(
      profile,
      profile[1L, ]
    )
    
    testthat::expect_error(
      mis_sek(
        profile = profile,
        global_reject = TRUE,
        c_n = 0.01,
        support_sets = list(
          overlap = 1:3,
          local = 1:3
        )
      ),
      "Duplicated pairs"
    )
  }
)


testthat::test_that(
  "factor profile columns are rejected",
  {
    profile <- make_valid_profile()
    profile$k <- factor(
      profile$k
    )
    
    testthat::expect_error(
      mis_sek(
        profile = profile,
        global_reject = TRUE,
        c_n = 0.01,
        support_sets = list(
          overlap = 1:3,
          local = 1:3
        )
      ),
      "silent factor"
    )
  }
)

# ------------------------------------------------------------------------------
# 3. Denominator admissibility and primary support
# ------------------------------------------------------------------------------

testthat::test_that(
  "denominator screening occurs before K0 construction",
  {
    profile <- make_valid_profile(
      k_grid = 1:2,
      plus = c(10, 1),
      minus = c(9, 0)
    )
    
    profile$denominator_fraction[
      profile$k == 1
    ] <- 0.01
    
    result <- mis_sek(
      profile = profile,
      global_reject = TRUE,
      c_n = 0.01,
      support_sets = list(
        overlap = 1:2,
        local = 1:2
      ),
      minimum_denominator_fraction = 0.05,
      denominator_radius = 0.01
    )
    
    testthat::expect_false(
      1L %in%
        result$primary$support_set
    )
    
    testthat::expect_true(
      2L %in%
        result$primary$support_set
    )
  }
)

# ------------------------------------------------------------------------------
# 4. Direction certification and set-valued selection
# ------------------------------------------------------------------------------

testthat::test_that(
  "direction requires a statistical gap larger than 2 c_n",
  {
    profile <- make_valid_profile(
      k_grid = 1,
      plus = 1.10,
      minus = 1.00
    )
    
    result <- mis_sek(
      profile = profile,
      global_reject = TRUE,
      c_n = 0.10,
      support_sets = list(
        overlap = 1,
        local = 1
      ),
      denominator_radius = 0.01
    )
    
    testthat::expect_identical(
      result$direction$support_set,
      integer(0L)
    )
  }
)


testthat::test_that(
  "set-valued output retains the common direction",
  {
    profile <- make_valid_profile(
      k_grid = 1:2,
      plus = c(1, 1),
      minus = c(0, 0)
    )
    
    result <- mis_sek(
      profile = profile,
      global_reject = TRUE,
      c_n = 0.01,
      support_sets = list(
        overlap = 1:2,
        local = 1:2
      ),
      denominator_radius = 0.01
    )
    
    testthat::expect_identical(
      result$selection_type,
      "set"
    )
    
    testthat::expect_identical(
      result$selected_sgn,
      1L
    )
  }
)

testthat::test_that(
  "row permutation cannot change the result",
  {
    profile <- make_valid_profile()
    
    result_a <- mis_sek(
      profile = profile,
      global_reject = TRUE,
      c_n = 0.01,
      support_sets = list(
        overlap = 1:3,
        local = 1:3
      ),
      denominator_radius = 0.01
    )
    
    set.seed(1)
    permuted <- profile[
      sample(
        seq_len(
          nrow(profile)
        )
      ),
      ,
      drop = FALSE
    ]
    
    result_b <- mis_sek(
      profile = permuted,
      global_reject = TRUE,
      c_n = 0.01,
      support_sets = list(
        overlap = 1:3,
        local = 1:3
      ),
      denominator_radius = 0.01
    )
    
    testthat::expect_identical(
      result_a$identified$k_set,
      result_b$identified$k_set
    )
    
    testthat::expect_identical(
      result_a$selected_sgn,
      result_b$selected_sgn
    )
  }
)

# ------------------------------------------------------------------------------
# 5. Determinism and provenance
# ------------------------------------------------------------------------------

testthat::test_that(
  "row permutation cannot change the result",
  {
    profile <- make_valid_profile()
    
    result_a <- mis_sek(
      profile = profile,
      global_reject = TRUE,
      c_n = 0.01,
      support_sets = list(
        overlap = 1:3,
        local = 1:3
      ),
      denominator_radius = 0.01
    )
    
    set.seed(1)
    permuted <- profile[
      sample(
        seq_len(
          nrow(profile)
        )
      ),
      ,
      drop = FALSE
    ]
    
    result_b <- mis_sek(
      profile = permuted,
      global_reject = TRUE,
      c_n = 0.01,
      support_sets = list(
        overlap = 1:3,
        local = 1:3
      ),
      denominator_radius = 0.01
    )
    
    testthat::expect_identical(
      result_a$identified$k_set,
      result_b$identified$k_set
    )
    
    testthat::expect_identical(
      result_a$selected_sgn,
      result_b$selected_sgn
    )
  }
)

testthat::test_that(
  "trace key comes from the exact selected sign row",
  {
    profile <- make_valid_profile(
      k_grid = 1,
      plus = 2,
      minus = 0
    )
    
    profile$trace_key[
      profile$sgn == 1
    ] <- "plus-row"
    
    profile$trace_key[
      profile$sgn == -1
    ] <- "minus-row"
    
    result <- mis_sek(
      profile = profile,
      global_reject = TRUE,
      c_n = 0.01,
      support_sets = list(
        overlap = 1,
        local = 1
      ),
      denominator_radius = 0.01,
      allow_automatic_submission = TRUE
    )
    
    testthat::expect_identical(
      result$selected_trace_key,
      "plus-row"
    )
    
    testthat::expect_true(
      result$automatic_top_k_submission
    )
  }
)

# ------------------------------------------------------------------------------
# 6. Set-intersection monotonicity
# ------------------------------------------------------------------------------

testthat::test_that(
  "tightening a support set cannot enlarge the final intersection",
  {
    profile <- make_valid_profile()
    
    broad <- mis_sek(
      profile = profile,
      global_reject = TRUE,
      c_n = 0.01,
      support_sets = list(
        overlap = 1:3,
        local = 1:3
      ),
      denominator_radius = 0.01
    )
    
    narrow <- mis_sek(
      profile = profile,
      global_reject = TRUE,
      c_n = 0.01,
      support_sets = list(
        overlap = 1:2,
        local = 1:3
      ),
      denominator_radius = 0.01
    )
    
    testthat::expect_true(
      all(
        narrow$identified$k_set %in%
          broad$identified$k_set
      )
    )
  }
)