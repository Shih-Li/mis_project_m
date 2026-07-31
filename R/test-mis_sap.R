# ==============================================================================
# File: /R/test-sap.R
# Purpose: Regression, mathematical-identity, calibration, and public-contract
#          tests for MIS-SAP influential-coalition analysis.
#
# Verifies:
#   - strict input and target validation;
#   - bounded adaptive-grid construction;
#   - exact deleted-refit coefficient identities;
#   - complete-search permutation calibration;
#   - deterministic same-seed behavior;
#   - sensitivity/cleaning separation;
#   - robust-anchor selection and safe abstention;
#   - location, scale, sign, and row-order properties;
#   - structured S3 output and print behavior;
#   - optional parity with the validated v2 wrapper during migration.
#
# Dependencies:
#   - /R/helpers_local.R
#   - /R/dinkelbach_topk.R
#   - /R/mis_sap.R
#   - testthat
#
# Optional migration dependency:
#   - /R/iterative_peel_sap.R
#   - /R/sap_exact_refit_v2.R
#   - /R/iterative_peel_sap_v2.R
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Test fixtures and comparison helpers
# ------------------------------------------------------------------------------

make_clean_sap_data <- function(
    n = 100L
) {
  x <- seq(
    -2,
    2,
    length.out = n
  )
  
  z <- cos(
    seq(
      0,
      3 * pi,
      length.out = n
    )
  )
  
  y <-
    1.25 +
    1.75 *
    x +
    0.40 *
    z +
    0.04 *
    sin(
      seq(
        0,
        5 * pi,
        length.out = n
      )
    )
  
  data.frame(
    y = y,
    x = x,
    z = z
  )
}


make_noiseless_sap_data <- function(
    n = 100L
) {
  x <- seq(
    -3,
    3,
    length.out = n
  )
  
  data.frame(
    y =
      2 +
      1.5 *
      x,
    x = x
  )
}


make_bad_leverage_sap_data <- function(
    n = 120L,
    k = 6L,
    seed = 731001L
) {
  stopifnot(
    n >= 40L,
    k >= 1L,
    k <= floor(
      0.05 *
        n
    )
  )
  
  set.seed(
    seed
  )
  
  z <- stats::rnorm(n)
  x <- stats::rnorm(n)
  
  y <-
    0.75 +
    1.80 *
    x +
    0.45 *
    z +
    stats::rnorm(
      n,
      sd = 0.35
    )
  
  truth <- seq_len(k)
  
  # Deterministic, same-direction, severe bad-leverage coalition.
  x[
    truth
  ] <-
    12 +
    seq_len(k) /
    100
  
  y[
    truth
  ] <-
    -85 -
    seq_len(k) /
    10
  
  result <- data.frame(
    y = y,
    x = x,
    z = z
  )
  
  attr(
    result,
    "truth"
  ) <- as.integer(
    truth
  )
  
  result
}


run_mis_sap_test <- function(
    data,
    seed = 731101L,
    target = "x",
    target_pos = NULL,
    k_grid = c(
      1L,
      2L,
      5L,
      6L
    ),
    B_perm = 19L,
    alpha = 0.05,
    max_fraction = 0.05,
    use_robust_anchor = FALSE,
    anchor_coefficient = NULL
) {
  set.seed(
    seed
  )
  
  mis_sap(
    formula = y ~ x + z,
    data = data,
    target = target,
    target_pos = target_pos,
    k_grid = k_grid,
    B_perm = B_perm,
    alpha = alpha,
    max_fraction = max_fraction,
    use_robust_anchor =
      use_robust_anchor,
    anchor_coefficient =
      anchor_coefficient
  )
}


expect_detection_equal <- function(
    first,
    second,
    tolerance = 1e-12
) {
  testthat::expect_identical(
    first$global$reject,
    second$global$reject
  )
  
  testthat::expect_equal(
    first$global$p_value,
    second$global$p_value,
    tolerance = tolerance
  )
  
  testthat::expect_equal(
    first$global$statistic,
    second$global$statistic,
    tolerance = tolerance
  )
  
  testthat::expect_identical(
    first$sensitivity$selected_k,
    second$sensitivity$selected_k
  )
  
  testthat::expect_identical(
    first$sensitivity$direction,
    second$sensitivity$direction
  )
  
  testthat::expect_identical(
    sort(
      first$sensitivity$set
    ),
    sort(
      second$sensitivity$set
    )
  )
  
  testthat::expect_equal(
    first$profile,
    second$profile,
    tolerance = tolerance
  )
}


# ------------------------------------------------------------------------------
# 1. Public interface, state, and structured-result contract
# ------------------------------------------------------------------------------

testthat::test_that(
  "mis_sap returns the formal structured S3 object",
  {
    data_case <- make_clean_sap_data()
    
    set.seed(731201L)
    
    result <- mis_sap(
      formula = y ~ x + z,
      data = data_case,
      target = "x",
      k_grid = c(
        1L,
        2L,
        5L
      ),
      B_perm = 19L,
      use_robust_anchor = FALSE
    )
    
    testthat::expect_s3_class(
      result,
      "mis_sap"
    )
    
    testthat::expect_s3_class(
      result,
      "list"
    )
    
    testthat::expect_identical(
      result$method,
      "MIS-SAP"
    )
    
    testthat::expect_identical(
      result$version,
      "2.0.0"
    )
    
    testthat::expect_true(
      result$state %in%
        c(
          "no_excessive_influence",
          "sensitivity_only",
          "robust_cleaning_available"
        )
    )
    
    testthat::expect_identical(
      result$target$name,
      "x"
    )
    
    testthat::expect_identical(
      result$target$position,
      2L
    )
    
    testthat::expect_named(
      result,
      c(
        "method",
        "version",
        "state",
        "target",
        "global",
        "grid",
        "sensitivity",
        "cleaning",
        "directional_candidates",
        "profile",
        "assumptions",
        "call"
      ),
      ignore.order = FALSE
    )
    
    testthat::expect_named(
      result$global,
      c(
        "reject",
        "p_value",
        "alpha",
        "statistic",
        "null_q95",
        "peak_excess_ratio",
        "permutations",
        "minimum_attainable_p",
        "stop_reason"
      ),
      ignore.order = FALSE
    )
    
    testthat::expect_named(
      result$sensitivity,
      c(
        "detected",
        "selected_k",
        "direction",
        "set",
        "data_rows",
        "row_names",
        "exact"
      ),
      ignore.order = FALSE
    )
    
    testthat::expect_named(
      result$cleaning,
      c(
        "permitted",
        "direction",
        "set",
        "data_rows",
        "row_names",
        "coefficient_retained",
        "anchor_coefficient",
        "anchor_source",
        "anchor_valid",
        "candidate_valid",
        "tied",
        "status",
        "positive_set",
        "negative_set",
        "positive_coefficient",
        "negative_coefficient",
        "positive_distance",
        "negative_distance",
        "anchor_error",
        "positive_error",
        "negative_error"
      ),
      ignore.order = FALSE
    )
    
    testthat::expect_identical(
      result$global$permutations,
      19L
    )
    
    testthat::expect_equal(
      result$global$minimum_attainable_p,
      1 / 20,
      tolerance = 0
    )
    
    testthat::expect_identical(
      result$grid$source,
      "explicit"
    )
    
    testthat::expect_identical(
      result$grid$values,
      c(
        1L,
        2L,
        5L
      )
    )
  }
)


testthat::test_that(
  "target name and target position resolve to the same analysis",
  {
    data_case <- make_clean_sap_data()
    
    set.seed(731202L)
    
    by_name <- mis_sap(
      y ~ x + z,
      data = data_case,
      target = "x",
      k_grid = c(
        1L,
        2L,
        5L
      ),
      B_perm = 19L,
      use_robust_anchor = FALSE
    )
    
    set.seed(731202L)
    
    by_position <- mis_sap(
      y ~ x + z,
      data = data_case,
      target_pos = 2L,
      k_grid = c(
        1L,
        2L,
        5L
      ),
      B_perm = 19L,
      use_robust_anchor = FALSE
    )
    
    expect_detection_equal(
      by_name,
      by_position
    )
    
    testthat::expect_identical(
      by_name$target,
      by_position$target
    )
  }
)


testthat::test_that(
  "model-row and original-data row positions are explicitly retained",
  {
    data_case <- make_bad_leverage_sap_data()
    
    row.names(data_case) <- paste0(
      "case-",
      seq_len(
        nrow(data_case)
      )
    )
    
    result <- run_mis_sap_test(
      data = data_case,
      seed = 731203L
    )
    
    testthat::expect_true(
      result$global$reject
    )
    
    testthat::expect_identical(
      result$sensitivity$data_rows,
      result$sensitivity$set
    )
    
    testthat::expect_identical(
      result$sensitivity$row_names,
      row.names(data_case)[
        result$sensitivity$set
      ]
    )
  }
)


# ------------------------------------------------------------------------------
# 2. Strict public-input and control validation
# ------------------------------------------------------------------------------

testthat::test_that(
  "a target must be supplied and must identify one non-intercept coefficient",
  {
    data_case <- make_clean_sap_data()
    
    testthat::expect_error(
      mis_sap(
        y ~ x + z,
        data = data_case,
        k_grid = 1L,
        B_perm = 19L
      ),
      "Supply either"
    )
    
    testthat::expect_error(
      mis_sap(
        y ~ x + z,
        data = data_case,
        target = "missing",
        k_grid = 1L,
        B_perm = 19L
      ),
      "must exactly identify"
    )
    
    testthat::expect_error(
      mis_sap(
        y ~ x + z,
        data = data_case,
        target = "(Intercept)",
        k_grid = 1L,
        B_perm = 19L
      ),
      "intercept"
    )
    
    testthat::expect_error(
      mis_sap(
        y ~ x + z,
        data = data_case,
        target = "x",
        target_pos = 3L,
        k_grid = 1L,
        B_perm = 19L
      ),
      "different coefficients"
    )
  }
)


testthat::test_that(
  "permutation and analysis controls are strictly validated",
  {
    testthat::expect_error(
      .sap_validate_test_controls(
        B_perm = 18L,
        alpha = 0.05,
        max_fraction = 0.05,
        use_robust_anchor = TRUE,
        tie_tolerance = 1e-10
      ),
      "too small"
    )
    
    testthat::expect_error(
      .sap_validate_test_controls(
        B_perm = 19.5,
        alpha = 0.05,
        max_fraction = 0.05,
        use_robust_anchor = TRUE,
        tie_tolerance = 1e-10
      ),
      "`B_perm`"
    )
    
    testthat::expect_error(
      .sap_validate_test_controls(
        B_perm = 19L,
        alpha = 0,
        max_fraction = 0.05,
        use_robust_anchor = TRUE,
        tie_tolerance = 1e-10
      ),
      "`alpha`"
    )
    
    testthat::expect_error(
      .sap_validate_test_controls(
        B_perm = 19L,
        alpha = 0.05,
        max_fraction = 0,
        use_robust_anchor = TRUE,
        tie_tolerance = 1e-10
      ),
      "`max_fraction`"
    )
    
    testthat::expect_error(
      .sap_validate_test_controls(
        B_perm = 19L,
        alpha = 0.05,
        max_fraction = 0.05,
        use_robust_anchor = NA,
        tie_tolerance = 1e-10
      ),
      "`use_robust_anchor`"
    )
    
    testthat::expect_error(
      .sap_validate_test_controls(
        B_perm = 19L,
        alpha = 0.05,
        max_fraction = 0.05,
        use_robust_anchor = TRUE,
        tie_tolerance = -1
      ),
      "`tie_tolerance`"
    )
    
    testthat::expect_silent(
      .sap_validate_test_controls(
        B_perm = 19L,
        alpha = 0.05,
        max_fraction = 0.05,
        use_robust_anchor = FALSE,
        tie_tolerance = 0
      )
    )
  }
)


testthat::test_that(
  "nonfinite data and inadmissible explicit grids are rejected",
  {
    data_case <- make_clean_sap_data()
    
    data_bad <- data_case
    data_bad$y[1L] <- NA_real_
    
    testthat::expect_error(
      mis_sap(
        y ~ x + z,
        data = data_bad,
        target = "x",
        k_grid = 1L,
        B_perm = 19L
      )
    )
    
    testthat::expect_error(
      mis_sap(
        y ~ x + z,
        data = data_case,
        target = "x",
        k_grid = c(
          20L,
          30L
        ),
        B_perm = 19L,
        max_fraction = 0.05
      ),
      "No explicit candidate"
    )
    
    testthat::expect_error(
      .sap_validate_delete_indices(
        delete = c(
          0L,
          1L
        ),
        n = 10L
      ),
      "valid integer row positions"
    )
  }
)


# ------------------------------------------------------------------------------
# 3. Bounded adaptive-grid contract
# ------------------------------------------------------------------------------

testthat::test_that(
  "adaptive grids reproduce the validated representative grids",
  {
    testthat::expect_identical(
      .sap_adaptive_k_grid(
        n = 500L,
        design_columns = 2L
      ),
      c(
        1L,
        2L,
        5L,
        10L,
        12L,
        20L,
        25L
      )
    )
    
    testthat::expect_identical(
      .sap_adaptive_k_grid(
        n = 1000L,
        design_columns = 2L
      ),
      c(
        1L,
        2L,
        5L,
        10L,
        20L,
        25L,
        50L
      )
    )
    
    testthat::expect_identical(
      .sap_adaptive_k_grid(
        n = 2500L,
        design_columns = 2L
      ),
      c(
        1L,
        2L,
        5L,
        10L,
        12L,
        20L,
        25L,
        62L,
        125L
      )
    )
    
    testthat::expect_identical(
      .sap_adaptive_k_grid(
        n = 5000L,
        design_columns = 2L
      ),
      c(
        1L,
        2L,
        5L,
        10L,
        20L,
        25L,
        50L,
        125L,
        250L
      )
    )
  }
)


testthat::test_that(
  "adaptive grids are sorted, unique, positive, and correctly capped",
  {
    sample_sizes <- c(
      20L,
      seq(
        25L,
        5000L,
        by = 25L
      ),
      10000L,
      25000L
    )
    
    for (n in sample_sizes) {
      grid <- .sap_adaptive_k_grid(
        n = n,
        design_columns = 2L
      )
      
      testthat::expect_identical(
        grid,
        sort(
          unique(grid)
        )
      )
      
      testthat::expect_true(
        all(
          grid >= 1L
        )
      )
      
      testthat::expect_true(
        all(
          grid <=
            floor(
              0.05 *
                n
            )
        )
      )
      
      testthat::expect_true(
        all(
          grid <=
            n -
            2L -
            1L
        )
      )
    }
  }
)


testthat::test_that(
  "explicit grids are normalized without exceeding the validated cap",
  {
    result <- .sap_validate_explicit_k_grid(
      k_grid = c(
        5L,
        1L,
        5L,
        10L
      ),
      n = 100L,
      design_columns = 2L,
      max_fraction = 0.05
    )
    
    testthat::expect_identical(
      result,
      c(
        1L,
        5L
      )
    )
  }
)


# ------------------------------------------------------------------------------
# 4. Exact deleted-refit identities
# ------------------------------------------------------------------------------

testthat::test_that(
  "exact simple-regression deletion equals direct OLS refitting",
  {
    data_case <- make_clean_sap_data()
    
    components <- .sap_build_model_components(
      y ~ x,
      data = data_case
    )
    
    delete <- c(
      1L,
      4L,
      9L,
      17L
    )
    
    exact <- .sap_exact_refit_target(
      components = components,
      target_pos = 2L,
      delete = delete
    )
    
    direct_full <- unname(
      stats::coef(
        stats::lm(
          y ~ x,
          data = data_case
        )
      )[2L]
    )
    
    direct_retained <- unname(
      stats::coef(
        stats::lm(
          y ~ x,
          data = data_case[
            -delete,
            ,
            drop = FALSE
          ]
        )
      )[2L]
    )
    
    testthat::expect_equal(
      exact$coefficient_full,
      direct_full,
      tolerance = 1e-10
    )
    
    testthat::expect_equal(
      exact$coefficient_retained,
      direct_retained,
      tolerance = 1e-10
    )
    
    testthat::expect_equal(
      exact$exact_shift,
      direct_retained -
        direct_full,
      tolerance = 1e-10
    )
    
    testthat::expect_true(
      exact$denominator_fraction >
        0
    )
  }
)


testthat::test_that(
  "exact nuisance-adjusted deletion equals direct multivariable OLS",
  {
    set.seed(731401L)
    
    n <- 160L
    z1 <- stats::rnorm(n)
    z2 <- stats::rnorm(n)
    
    x <-
      0.4 *
      z1 -
      0.3 *
      z2 +
      stats::rnorm(n)
    
    y <-
      0.5 +
      1.9 *
      x +
      0.6 *
      z1 -
      0.2 *
      z2 +
      stats::rnorm(
        n,
        sd = 0.7
      )
    
    data_case <- data.frame(
      y = y,
      x = x,
      z1 = z1,
      z2 = z2
    )
    
    delete <- c(
      2L,
      5L,
      11L,
      40L,
      91L
    )
    
    components <- .sap_build_model_components(
      y ~ x + z1 + z2,
      data = data_case
    )
    
    exact <- .sap_exact_refit_target(
      components = components,
      target_pos = 2L,
      delete = delete
    )
    
    direct_full <- unname(
      stats::coef(
        stats::lm(
          y ~ x + z1 + z2,
          data = data_case
        )
      )[2L]
    )
    
    direct_retained <- unname(
      stats::coef(
        stats::lm(
          y ~ x + z1 + z2,
          data = data_case[
            -delete,
            ,
            drop = FALSE
          ]
        )
      )[2L]
    )
    
    testthat::expect_equal(
      exact$coefficient_full,
      direct_full,
      tolerance = 1e-10
    )
    
    testthat::expect_equal(
      exact$coefficient_retained,
      direct_retained,
      tolerance = 1e-10
    )
    
    testthat::expect_equal(
      exact$exact_shift,
      direct_retained -
        direct_full,
      tolerance = 1e-10
    )
  }
)


testthat::test_that(
  "exact refitting is invariant to row ordering",
  {
    set.seed(731402L)
    
    n <- 120L
    z <- stats::rnorm(n)
    x <- stats::rnorm(n)
    
    y <-
      0.8 +
      1.7 *
      x +
      0.5 *
      z +
      stats::rnorm(n)
    
    data_original <- data.frame(
      y = y,
      x = x,
      z = z
    )
    
    delete_original <- c(
      3L,
      14L,
      48L,
      79L
    )
    
    original <- .sap_exact_refit_target(
      components =
        .sap_build_model_components(
          y ~ x + z,
          data = data_original
        ),
      target_pos = 2L,
      delete = delete_original
    )
    
    permutation <- sample.int(n)
    
    data_permuted <- data_original[
      permutation,
      ,
      drop = FALSE
    ]
    
    delete_permuted <- match(
      delete_original,
      permutation
    )
    
    permuted <- .sap_exact_refit_target(
      components =
        .sap_build_model_components(
          y ~ x + z,
          data = data_permuted
        ),
      target_pos = 2L,
      delete = delete_permuted
    )
    
    testthat::expect_equal(
      original$coefficient_full,
      permuted$coefficient_full,
      tolerance = 1e-10
    )
    
    testthat::expect_equal(
      original$coefficient_retained,
      permuted$coefficient_retained,
      tolerance = 1e-10
    )
    
    testthat::expect_equal(
      original$exact_shift,
      permuted$exact_shift,
      tolerance = 1e-10
    )
  }
)


testthat::test_that(
  "degenerate retained target energy is rejected",
  {
    data_case <- data.frame(
      x = c(
        rep(
          0,
          18L
        ),
        1,
        2
      )
    )
    
    data_case$y <-
      1 +
      2 *
      data_case$x
    
    components <- .sap_build_model_components(
      y ~ x,
      data = data_case
    )
    
    testthat::expect_error(
      .sap_exact_refit_target(
        components = components,
        target_pos = 2L,
        delete = c(
          19L,
          20L
        )
      ),
      "insufficient residualized energy"
    )
  }
)


# ------------------------------------------------------------------------------
# 5. Deterministic directional search and row-order invariance
# ------------------------------------------------------------------------------

testthat::test_that(
  "directional search retains both candidates and is deterministic",
  {
    data_case <- make_clean_sap_data(
      n = 120L
    )
    
    components <- .sap_build_model_components(
      y ~ x + z,
      data = data_case
    )
    
    fitted <- .sap_fwl_components(
      X = components$X,
      y = components$y,
      target_pos = 2L
    )
    
    first <- .sap_directional_candidates_at_k(
      fwl_fit = fitted,
      k = 6L
    )
    
    second <- .sap_directional_candidates_at_k(
      fwl_fit = fitted,
      k = 6L
    )
    
    testthat::expect_identical(
      first,
      second
    )
    
    testthat::expect_length(
      first$positive_set,
      6L
    )
    
    testthat::expect_length(
      first$negative_set,
      6L
    )
    
    testthat::expect_length(
      first$sensitivity_set,
      6L
    )
    
    testthat::expect_true(
      first$sensitivity_direction %in%
        c(
          -1L,
          1L
        )
    )
  }
)


testthat::test_that(
  "directional coalition search is invariant to row ordering",
  {
    set.seed(731501L)
    
    n <- 120L
    z <- stats::rnorm(n)
    
    x <-
      0.4 *
      z +
      stats::rnorm(n)
    
    y <-
      0.5 +
      1.6 *
      x +
      0.5 *
      z +
      stats::rnorm(n)
    
    data_original <- data.frame(
      y = y,
      x = x,
      z = z
    )
    
    original_components <- .sap_build_model_components(
      y ~ x + z,
      data = data_original
    )
    
    original_fit <- .sap_fwl_components(
      X = original_components$X,
      y = original_components$y,
      target_pos = 2L
    )
    
    original <- .sap_directional_candidates_at_k(
      fwl_fit = original_fit,
      k = 6L
    )
    
    permutation <- sample.int(n)
    
    data_permuted <- data_original[
      permutation,
      ,
      drop = FALSE
    ]
    
    permuted_components <- .sap_build_model_components(
      y ~ x + z,
      data = data_permuted
    )
    
    permuted_fit <- .sap_fwl_components(
      X = permuted_components$X,
      y = permuted_components$y,
      target_pos = 2L
    )
    
    permuted <- .sap_directional_candidates_at_k(
      fwl_fit = permuted_fit,
      k = 6L
    )
    
    testthat::expect_identical(
      sort(
        original$positive_set
      ),
      sort(
        permutation[
          permuted$positive_set
        ]
      )
    )
    
    testthat::expect_identical(
      sort(
        original$negative_set
      ),
      sort(
        permutation[
          permuted$negative_set
        ]
      )
    )
    
    testthat::expect_identical(
      sort(
        original$sensitivity_set
      ),
      sort(
        permutation[
          permuted$sensitivity_set
        ]
      )
    )
    
    testthat::expect_identical(
      original$sensitivity_direction,
      permuted$sensitivity_direction
    )
    
    testthat::expect_equal(
      original$positive_search_shift,
      permuted$positive_search_shift,
      tolerance = 1e-10
    )
    
    testthat::expect_equal(
      original$negative_search_shift,
      permuted$negative_search_shift,
      tolerance = 1e-10
    )
  }
)


# ------------------------------------------------------------------------------
# 6. Non-rejection safety and reproducibility
# ------------------------------------------------------------------------------

testthat::test_that(
  "a noiseless exact model creates no sensitivity or cleaning deletion",
  {
    data_case <- make_noiseless_sap_data()
    
    set.seed(2160001L)
    
    result <- mis_sap(
      y ~ x,
      data = data_case,
      target = "x",
      k_grid = c(
        1L,
        2L,
        5L
      ),
      B_perm = 39L,
      alpha = 0.05,
      use_robust_anchor = TRUE
    )
    
    testthat::expect_false(
      result$global$reject
    )
    
    testthat::expect_identical(
      result$state,
      "no_excessive_influence"
    )
    
    testthat::expect_identical(
      result$sensitivity$selected_k,
      0L
    )
    
    testthat::expect_identical(
      result$sensitivity$set,
      integer(0L)
    )
    
    testthat::expect_null(
      result$sensitivity$exact
    )
    
    testthat::expect_false(
      result$cleaning$permitted
    )
    
    testthat::expect_identical(
      result$cleaning$set,
      integer(0L)
    )
    
    testthat::expect_identical(
      result$cleaning$status,
      "not_detected_no_deletion"
    )
  }
)


testthat::test_that(
  "the complete result is exactly reproducible with the same seed",
  {
    data_case <- make_bad_leverage_sap_data()
    
    first <- run_mis_sap_test(
      data = data_case,
      seed = 731601L
    )
    
    second <- run_mis_sap_test(
      data = data_case,
      seed = 731601L
    )
    
    testthat::expect_identical(
      first,
      second
    )
  }
)


testthat::test_that(
  "permutation p-values obey the finite-sample counting grid",
  {
    data_case <- make_bad_leverage_sap_data()
    
    result <- run_mis_sap_test(
      data = data_case,
      seed = 731602L,
      B_perm = 19L
    )
    
    scaled_global_p <-
      result$global$p_value *
      (
        result$global$permutations +
          1L
      )
    
    scaled_profile_p <-
      result$profile$
      p_selection_adjusted *
      (
        result$global$permutations +
          1L
      )
    
    testthat::expect_equal(
      scaled_global_p,
      round(
        scaled_global_p
      ),
      tolerance = 1e-12
    )
    
    testthat::expect_equal(
      scaled_profile_p,
      round(
        scaled_profile_p
      ),
      tolerance = 1e-12
    )
    
    testthat::expect_true(
      result$global$p_value >=
        result$global$
        minimum_attainable_p
    )
  }
)


# ------------------------------------------------------------------------------
# 7. Exact post-detection sensitivity reporting
# ------------------------------------------------------------------------------

testthat::test_that(
  "a severe bad-leverage coalition is detected",
  {
    data_case <- make_bad_leverage_sap_data()
    
    result <- run_mis_sap_test(
      data = data_case,
      seed = 731701L
    )
    
    testthat::expect_true(
      result$global$reject
    )
    
    testthat::expect_identical(
      result$global$p_value,
      1 / 20
    )
    
    testthat::expect_true(
      result$sensitivity$selected_k %in%
        c(
          1L,
          2L,
          5L,
          6L
        )
    )
    
    testthat::expect_length(
      result$sensitivity$set,
      result$sensitivity$selected_k
    )
    
    testthat::expect_true(
      result$sensitivity$direction %in%
        c(
          -1L,
          1L
        )
    )
    
    testthat::expect_true(
      is.list(
        result$sensitivity$exact
      )
    )
  }
)


testthat::test_that(
  "reported exact sensitivity coefficient equals direct OLS refitting",
  {
    data_case <- make_bad_leverage_sap_data()
    
    result <- run_mis_sap_test(
      data = data_case,
      seed = 731702L
    )
    
    testthat::expect_true(
      result$global$reject
    )
    
    direct_full <- unname(
      stats::coef(
        stats::lm(
          y ~ x + z,
          data = data_case
        )
      )[2L]
    )
    
    direct_retained <- unname(
      stats::coef(
        stats::lm(
          y ~ x + z,
          data = data_case[
            -result$sensitivity$set,
            ,
            drop = FALSE
          ]
        )
      )[2L]
    )
    
    testthat::expect_equal(
      result$sensitivity$exact$
        coefficient_full,
      direct_full,
      tolerance = 1e-10
    )
    
    testthat::expect_equal(
      result$sensitivity$exact$
        coefficient_retained,
      direct_retained,
      tolerance = 1e-10
    )
    
    testthat::expect_equal(
      result$sensitivity$exact$
        exact_shift,
      direct_retained -
        direct_full,
      tolerance = 1e-10
    )
  }
)


# ------------------------------------------------------------------------------
# 8. Robust-anchor direction and abstention contract
# ------------------------------------------------------------------------------

testthat::test_that(
  "a supplied anchor selects the mathematically closer exact coefficient",
  {
    set.seed(731801L)
    
    n <- 120L
    x <- stats::rnorm(n)
    z <- stats::rnorm(n)
    
    y <-
      0.5 +
      1.8 *
      x +
      0.4 *
      z +
      stats::rnorm(
        n,
        sd = 0.7
      )
    
    data_case <- data.frame(
      y = y,
      x = x,
      z = z
    )
    
    components <- .sap_build_model_components(
      y ~ x + z,
      data = data_case
    )
    
    positive_set <- order(
      x,
      decreasing = TRUE
    )[1:6]
    
    negative_set <- order(
      x,
      decreasing = FALSE
    )[1:6]
    
    positive_coefficient <-
      .sap_exact_refit_target(
        components = components,
        target_pos = 2L,
        delete = positive_set
      )$coefficient_retained
    
    result <- .sap_robust_anchor_direction(
      components = components,
      target_pos = 2L,
      positive_set = positive_set,
      negative_set = negative_set,
      anchor_coefficient =
        positive_coefficient
    )
    
    testthat::expect_identical(
      result$selected_direction,
      1L
    )
    
    testthat::expect_identical(
      result$selected_set,
      sort(
        as.integer(
          positive_set
        )
      )
    )
    
    testthat::expect_true(
      result$automatic_deletion_permitted
    )
    
    testthat::expect_identical(
      result$decision_status,
      "positive_direction"
    )
  }
)


testthat::test_that(
  "anchor ties and invalid anchors safely prevent automatic deletion",
  {
    set.seed(731802L)
    
    n <- 100L
    x <- stats::rnorm(n)
    y <-
      1 +
      2 *
      x +
      stats::rnorm(n)
    
    data_case <- data.frame(
      y = y,
      x = x
    )
    
    components <- .sap_build_model_components(
      y ~ x,
      data = data_case
    )
    
    positive_set <- order(
      x,
      decreasing = TRUE
    )[1:5]
    
    negative_set <- order(
      x,
      decreasing = FALSE
    )[1:5]
    
    positive_coefficient <-
      .sap_exact_refit_target(
        components = components,
        target_pos = 2L,
        delete = positive_set
      )$coefficient_retained
    
    negative_coefficient <-
      .sap_exact_refit_target(
        components = components,
        target_pos = 2L,
        delete = negative_set
      )$coefficient_retained
    
    midpoint_anchor <- (
      positive_coefficient +
        negative_coefficient
    ) /
      2
    
    tied <- .sap_robust_anchor_direction(
      components = components,
      target_pos = 2L,
      positive_set = positive_set,
      negative_set = negative_set,
      anchor_coefficient =
        midpoint_anchor
    )
    
    invalid <- .sap_robust_anchor_direction(
      components = components,
      target_pos = 2L,
      positive_set = positive_set,
      negative_set = negative_set,
      anchor_coefficient = NA_real_
    )
    
    testthat::expect_identical(
      tied$selected_direction,
      0L
    )
    
    testthat::expect_false(
      tied$automatic_deletion_permitted
    )
    
    testthat::expect_identical(
      tied$decision_status,
      "anchor_tie_no_deletion"
    )
    
    testthat::expect_identical(
      invalid$selected_direction,
      0L
    )
    
    testthat::expect_false(
      invalid$automatic_deletion_permitted
    )
    
    testthat::expect_identical(
      invalid$decision_status,
      "invalid_anchor_no_deletion"
    )
  }
)


testthat::test_that(
  "invalid exact candidate refits prevent automatic deletion",
  {
    data_case <- data.frame(
      x = c(
        rep(
          0,
          18L
        ),
        1,
        2
      )
    )
    
    data_case$y <-
      1 +
      2 *
      data_case$x
    
    components <- .sap_build_model_components(
      y ~ x,
      data = data_case
    )
    
    result <- .sap_robust_anchor_direction(
      components = components,
      target_pos = 2L,
      positive_set = c(
        19L,
        20L
      ),
      negative_set = c(
        1L,
        2L
      ),
      anchor_coefficient = 2
    )
    
    testthat::expect_identical(
      result$selected_direction,
      0L
    )
    
    testthat::expect_false(
      result$automatic_deletion_permitted
    )
    
    testthat::expect_identical(
      result$decision_status,
      "invalid_candidate_no_deletion"
    )
  }
)


testthat::test_that(
  "automatic lmrob anchoring returns a finite anchor on regular data",
  {
    testthat::skip_if_not_installed(
      "robustbase"
    )
    
    set.seed(731803L)
    
    n <- 150L
    x <- stats::rnorm(n)
    
    y <-
      1 +
      2 *
      x +
      stats::rnorm(n)
    
    data_case <- data.frame(
      y = y,
      x = x
    )
    
    components <- .sap_build_model_components(
      y ~ x,
      data = data_case
    )
    
    result <- .sap_robust_anchor_direction(
      components = components,
      target_pos = 2L,
      positive_set =
        order(
          x,
          decreasing = TRUE
        )[1:5],
      negative_set =
        order(
          x,
          decreasing = FALSE
        )[1:5],
      anchor_coefficient = NULL
    )
    
    testthat::expect_true(
      result$anchor_valid
    )
    
    testthat::expect_true(
      is.finite(
        result$anchor_coefficient
      )
    )
    
    testthat::expect_true(
      result$selected_direction %in%
        c(
          -1L,
          0L,
          1L
        )
    )
  }
)


# ------------------------------------------------------------------------------
# 9. Detection and robust cleaning are strictly separated
# ------------------------------------------------------------------------------

testthat::test_that(
  "enabling the robust anchor cannot change calibrated detection",
  {
    data_case <- make_bad_leverage_sap_data()
    
    without_anchor <- run_mis_sap_test(
      data = data_case,
      seed = 731901L,
      use_robust_anchor = FALSE
    )
    
    with_anchor <- run_mis_sap_test(
      data = data_case,
      seed = 731901L,
      use_robust_anchor = TRUE,
      anchor_coefficient = 1.8
    )
    
    expect_detection_equal(
      without_anchor,
      with_anchor
    )
    
    testthat::expect_identical(
      without_anchor$grid,
      with_anchor$grid
    )
    
    if (with_anchor$cleaning$permitted) {
      testthat::expect_true(
        identical(
          sort(
            with_anchor$cleaning$set
          ),
          sort(
            with_anchor$
              directional_candidates$
              positive_set
          )
        ) ||
          identical(
            sort(
              with_anchor$cleaning$set
            ),
            sort(
              with_anchor$
                directional_candidates$
                negative_set
            )
          )
      )
    } else {
      testthat::expect_identical(
        with_anchor$cleaning$set,
        integer(0L)
      )
    }
  }
)


testthat::test_that(
  "disabling robust anchoring preserves directional candidates but abstains",
  {
    data_case <- make_bad_leverage_sap_data()
    
    result <- run_mis_sap_test(
      data = data_case,
      seed = 731902L,
      use_robust_anchor = FALSE
    )
    
    testthat::expect_true(
      result$global$reject
    )
    
    testthat::expect_false(
      result$cleaning$permitted
    )
    
    testthat::expect_identical(
      result$cleaning$status,
      "robust_anchor_disabled"
    )
    
    testthat::expect_identical(
      result$cleaning$positive_set,
      result$directional_candidates$
        positive_set
    )
    
    testthat::expect_identical(
      result$cleaning$negative_set,
      result$directional_candidates$
        negative_set
    )
  }
)


# ------------------------------------------------------------------------------
# 10. Algebraic transformation properties
# ------------------------------------------------------------------------------

testthat::test_that(
  "detection obeys location, positive-scale, and sign transformations",
  {
    data_base <- make_bad_leverage_sap_data()
    
    base <- run_mis_sap_test(
      data = data_base,
      seed = 732001L
    )
    
    testthat::expect_true(
      base$global$reject
    )
    
    data_y_location <- transform(
      data_base,
      y = y + 100
    )
    
    data_x_location <- transform(
      data_base,
      x = x + 50
    )
    
    data_y_scale <- transform(
      data_base,
      y = y * 3.5
    )
    
    data_y_sign <- transform(
      data_base,
      y = -y
    )
    
    y_location <- run_mis_sap_test(
      data = data_y_location,
      seed = 732001L
    )
    
    x_location <- run_mis_sap_test(
      data = data_x_location,
      seed = 732001L
    )
    
    y_scale <- run_mis_sap_test(
      data = data_y_scale,
      seed = 732001L
    )
    
    y_sign <- run_mis_sap_test(
      data = data_y_sign,
      seed = 732001L
    )
    
    expect_detection_equal(
      base,
      y_location,
      tolerance = 1e-9
    )
    
    expect_detection_equal(
      base,
      x_location,
      tolerance = 1e-9
    )
    
    testthat::expect_identical(
      base$global$reject,
      y_scale$global$reject
    )
    
    testthat::expect_equal(
      base$global$p_value,
      y_scale$global$p_value,
      tolerance = 1e-12
    )
    
    testthat::expect_identical(
      base$sensitivity$selected_k,
      y_scale$sensitivity$selected_k
    )
    
    testthat::expect_identical(
      sort(
        base$sensitivity$set
      ),
      sort(
        y_scale$sensitivity$set
      )
    )
    
    testthat::expect_identical(
      base$sensitivity$direction,
      y_scale$sensitivity$direction
    )
    
    testthat::expect_equal(
      y_scale$global$statistic,
      3.5 *
        base$global$statistic,
      tolerance = 1e-9
    )
    
    testthat::expect_identical(
      base$global$reject,
      y_sign$global$reject
    )
    
    testthat::expect_equal(
      base$global$p_value,
      y_sign$global$p_value,
      tolerance = 1e-12
    )
    
    testthat::expect_identical(
      base$sensitivity$selected_k,
      y_sign$sensitivity$selected_k
    )
    
    testthat::expect_identical(
      sort(
        base$sensitivity$set
      ),
      sort(
        y_sign$sensitivity$set
      )
    )
    
    testthat::expect_identical(
      y_sign$sensitivity$direction,
      -base$sensitivity$direction
    )
    
    testthat::expect_equal(
      y_sign$global$statistic,
      base$global$statistic,
      tolerance = 1e-9
    )
  }
)


testthat::test_that(
  "exact deleted-refit coefficients obey the same transformations",
  {
    data_base <- make_bad_leverage_sap_data()
    
    components_base <- .sap_build_model_components(
      y ~ x + z,
      data = data_base
    )
    
    delete <- attr(
      data_base,
      "truth"
    )
    
    base <- .sap_exact_refit_target(
      components = components_base,
      target_pos = 2L,
      delete = delete
    )
    
    transformations <- list(
      y_location = list(
        data = transform(
          data_base,
          y = y + 100
        ),
        multiplier = 1
      ),
      x_location = list(
        data = transform(
          data_base,
          x = x + 50
        ),
        multiplier = 1
      ),
      y_scale = list(
        data = transform(
          data_base,
          y = y * 3.5
        ),
        multiplier = 3.5
      ),
      x_scale = list(
        data = transform(
          data_base,
          x = x * 2.5
        ),
        multiplier = 1 / 2.5
      ),
      y_sign = list(
        data = transform(
          data_base,
          y = -y
        ),
        multiplier = -1
      ),
      x_sign = list(
        data = transform(
          data_base,
          x = -x
        ),
        multiplier = -1
      )
    )
    
    for (
      transformation in transformations
    ) {
      transformed <- .sap_exact_refit_target(
        components =
          .sap_build_model_components(
            y ~ x + z,
            data =
              transformation$data
          ),
        target_pos = 2L,
        delete = delete
      )
      
      testthat::expect_equal(
        transformed$coefficient_retained,
        transformation$multiplier *
          base$coefficient_retained,
        tolerance = 1e-10
      )
      
      testthat::expect_equal(
        transformed$exact_shift,
        transformation$multiplier *
          base$exact_shift,
        tolerance = 1e-10
      )
      
      testthat::expect_equal(
        transformed$denominator_fraction,
        base$denominator_fraction,
        tolerance = 1e-10
      )
    }
  }
)


# ------------------------------------------------------------------------------
# 11. S3 print behavior
# ------------------------------------------------------------------------------

testthat::test_that(
  "print.mis_sap reports the formal state and returns invisibly",
  {
    data_case <- make_noiseless_sap_data()
    
    set.seed(2160001L)
    
    result <- mis_sap(
      y ~ x,
      data = data_case,
      target = "x",
      k_grid = c(
        1L,
        2L,
        5L
      ),
      B_perm = 39L,
      use_robust_anchor = TRUE
    )
    
    printed <- utils::capture.output(
      returned <- print(
        result
      )
    )
    
    testthat::expect_identical(
      returned,
      result
    )
    
    testthat::expect_true(
      any(
        grepl(
          "MIS-SAP influential-coalition analysis",
          printed,
          fixed = TRUE
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "State:",
          printed,
          fixed = TRUE
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "Target: x",
          printed,
          fixed = TRUE
        )
      )
    )
    
    testthat::expect_true(
      any(
        grepl(
          "no deletion is recommended",
          printed,
          fixed = TRUE
        )
      )
    )
  }
)


# ------------------------------------------------------------------------------
# 12. Optional migration parity with the validated v2 wrapper
# ------------------------------------------------------------------------------

testthat::test_that(
  "consolidated mis_sap detection matches the validated v2 wrapper",
  {
    if (
      !exists(
        "sap_multiscale_test_v2",
        mode = "function",
        inherits = TRUE
      )
    ) {
      testthat::skip(
        "The legacy v2 wrapper is not loaded; migration parity was skipped."
      )
    }
    
    data_case <- make_bad_leverage_sap_data()
    
    permutation_seed <- 732201L
    
    set.seed(
      permutation_seed
    )
    
    legacy <- sap_multiscale_test_v2(
      formula = y ~ x + z,
      data = data_case,
      target_pos = 2L,
      k_grid = c(
        1L,
        2L,
        5L,
        6L
      ),
      B_perm = 19L,
      alpha = 0.05,
      max_fraction = 0.05,
      use_robust_anchor = FALSE
    )
    
    set.seed(
      permutation_seed
    )
    
    consolidated <- mis_sap(
      formula = y ~ x + z,
      data = data_case,
      target = "x",
      k_grid = c(
        1L,
        2L,
        5L,
        6L
      ),
      B_perm = 19L,
      alpha = 0.05,
      max_fraction = 0.05,
      use_robust_anchor = FALSE
    )
    
    testthat::expect_identical(
      consolidated$global$reject,
      legacy$detected
    )
    
    testthat::expect_equal(
      consolidated$global$p_value,
      legacy$global_p,
      tolerance = 0
    )
    
    testthat::expect_equal(
      consolidated$global$statistic,
      legacy$observed_global_stat,
      tolerance = 0
    )
    
    testthat::expect_identical(
      consolidated$sensitivity$selected_k,
      legacy$selected_k
    )
    
    testthat::expect_identical(
      consolidated$sensitivity$direction,
      legacy$selected_direction
    )
    
    testthat::expect_identical(
      sort(
        consolidated$sensitivity$set
      ),
      sort(
        legacy$sensitivity_set
      )
    )
    
    testthat::expect_equal(
      consolidated$profile,
      legacy$profile,
      tolerance = 0
    )
  }
)