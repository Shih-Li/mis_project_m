
repo_root <- Sys.getenv(
  "MIS_SAP_REPO",
  unset = "/root/mis_project_m"
)

source(
  file.path(
    repo_root,
    "R",
    "sap_exact_refit_v2.R"
  ),
  local = TRUE
)

source(
  file.path(
    repo_root,
    "R",
    "iterative_peel_sap_v2.R"
  ),
  local = TRUE
)


testthat::test_that(
  "pure anchor rule selects the closer coefficient",
  {
    positive <- sap_choose_direction_from_anchor_v2(
      positive_coefficient = 1.1,
      negative_coefficient = 4.5,
      anchor_coefficient = 1
    )

    negative <- sap_choose_direction_from_anchor_v2(
      positive_coefficient = -3,
      negative_coefficient = 2.1,
      anchor_coefficient = 2
    )

    testthat::expect_identical(
      positive$selected_direction,
      1L
    )

    testthat::expect_identical(
      negative$selected_direction,
      -1L
    )
  }
)


testthat::test_that(
  "anchor ties lead to no automatic deletion",
  {
    result <- sap_choose_direction_from_anchor_v2(
      positive_coefficient = 1,
      negative_coefficient = 3,
      anchor_coefficient = 2
    )

    testthat::expect_identical(
      result$selected_direction,
      0L
    )

    testthat::expect_true(
      result$tied
    )

    testthat::expect_identical(
      result$decision_status,
      "anchor_tie_no_deletion"
    )
  }
)


testthat::test_that(
  "invalid anchors lead to no automatic deletion",
  {
    result <- sap_choose_direction_from_anchor_v2(
      positive_coefficient = 1,
      negative_coefficient = 3,
      anchor_coefficient = NA_real_
    )

    testthat::expect_identical(
      result$selected_direction,
      0L
    )

    testthat::expect_false(
      result$anchor_valid
    )

    testthat::expect_identical(
      result$decision_status,
      "invalid_anchor_no_deletion"
    )
  }
)


testthat::test_that(
  "high-level rule uses exact retained coefficients",
  {
    set.seed(44001)

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
        sd = 0.5
      )

    data_case <- data.frame(
      y = y,
      x = x,
      z = z
    )

    positive_set <- order(
      x,
      decreasing = TRUE
    )[1:6]

    negative_set <- order(
      x,
      decreasing = FALSE
    )[1:6]

    positive_direct <- stats::coef(
      stats::lm(
        y ~ x + z,
        data = data_case[
          -positive_set,
          ,
          drop = FALSE
        ]
      )
    )[2L]

    negative_direct <- stats::coef(
      stats::lm(
        y ~ x + z,
        data = data_case[
          -negative_set,
          ,
          drop = FALSE
        ]
      )
    )[2L]

    anchor <- positive_direct + 1e-4

    result <- sap_robust_anchor_direction_v2(
      formula = y ~ x + z,
      data = data_case,
      target_pos = 2L,
      positive_set = positive_set,
      negative_set = negative_set,
      anchor_coefficient = anchor
    )

    testthat::expect_equal(
      result$positive_coefficient,
      unname(
        positive_direct
      ),
      tolerance = 1e-10
    )

    testthat::expect_equal(
      result$negative_coefficient,
      unname(
        negative_direct
      ),
      tolerance = 1e-10
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
  }
)


testthat::test_that(
  "automatic lmrob anchors are finite on regular data",
  {
    set.seed(44002)

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

    positive_set <- order(
      x,
      decreasing = TRUE
    )[1:5]

    negative_set <- order(
      x,
      decreasing = FALSE
    )[1:5]

    result <- sap_robust_anchor_direction_v2(
      formula = y ~ x,
      data = data_case,
      target_pos = 2L,
      positive_set = positive_set,
      negative_set = negative_set
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


testthat::test_that(
  "row reordering preserves the selected original coalition",
  {
    set.seed(44003)

    n <- 100L

    z <- stats::rnorm(n)

    x <-
      0.3 *
        z +
      stats::rnorm(n)

    y <-
      0.8 +
      1.5 *
        x +
      0.5 *
        z +
      stats::rnorm(n)

    data_original <- data.frame(
      y = y,
      x = x,
      z = z
    )

    positive_original <- order(
      x,
      decreasing = TRUE
    )[1:5]

    negative_original <- order(
      x,
      decreasing = FALSE
    )[1:5]

    positive_coefficient <-
      sap_exact_refit_target_v2(
        y ~ x + z,
        data = data_original,
        target_pos = 2L,
        delete = positive_original
      )$coefficient_retained

    anchor <- positive_coefficient

    original_result <-
      sap_robust_anchor_direction_v2(
        formula = y ~ x + z,
        data = data_original,
        target_pos = 2L,
        positive_set = positive_original,
        negative_set = negative_original,
        anchor_coefficient = anchor
      )

    permutation <- sample.int(n)

    data_permuted <- data_original[
      permutation,
      ,
      drop = FALSE
    ]

    positive_permuted <- match(
      positive_original,
      permutation
    )

    negative_permuted <- match(
      negative_original,
      permutation
    )

    permuted_result <-
      sap_robust_anchor_direction_v2(
        formula = y ~ x + z,
        data = data_permuted,
        target_pos = 2L,
        positive_set = positive_permuted,
        negative_set = negative_permuted,
        anchor_coefficient = anchor
      )

    selected_permuted_original_rows <-
      sort(
        permutation[
          permuted_result$
            selected_set
        ]
      )

    testthat::expect_identical(
      original_result$selected_direction,
      permuted_result$selected_direction
    )

    testthat::expect_identical(
      sort(
        original_result$
          selected_set
      ),
      selected_permuted_original_rows
    )
  }
)


testthat::test_that(
  "invalid candidate refits prevent automatic deletion",
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

    result <- sap_robust_anchor_direction_v2(
      formula = y ~ x,
      data = data_case,
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

