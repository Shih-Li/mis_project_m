
repo_root <- Sys.getenv(
  "MIS_SAP_REPO",
  unset = "/root/mis_project_m"
)

for (
  path in file.path(
    repo_root,
    "R",
    c(
      "helpers_local.R",
      "dinkelbach_topk.R",
      "iterative_peel_sap.R",
      "sap_exact_refit_v2.R",
      "iterative_peel_sap_v2.R"
    )
  )
) {
  source(
    path,
    local = TRUE
  )
}


testthat::test_that(
  "permutation controls are validated",
  {
    testthat::expect_error(
      .sap_validate_test_controls_v2(
        B_perm = 18L,
        alpha = 0.05,
        use_robust_anchor = TRUE,
        tie_tolerance = 1e-10
      ),
      "too small"
    )

    testthat::expect_error(
      .sap_validate_test_controls_v2(
        B_perm = 19.5,
        alpha = 0.05,
        use_robust_anchor = TRUE,
        tie_tolerance = 1e-10
      ),
      "`B_perm`"
    )

    testthat::expect_error(
      .sap_validate_test_controls_v2(
        B_perm = 39L,
        alpha = 1,
        use_robust_anchor = TRUE,
        tie_tolerance = 1e-10
      ),
      "`alpha`"
    )

    testthat::expect_error(
      .sap_validate_test_controls_v2(
        B_perm = 39L,
        alpha = 0.05,
        use_robust_anchor = NA,
        tie_tolerance = 1e-10
      ),
      "`use_robust_anchor`"
    )

    testthat::expect_silent(
      .sap_validate_test_controls_v2(
        B_perm = 19L,
        alpha = 0.05,
        use_robust_anchor = FALSE,
        tie_tolerance = 0
      )
    )
  }
)


testthat::test_that(
  "same data and seed reproduce the complete v2 result",
  {
    set.seed(66001)

    n <- 120L

    x <- stats::rnorm(n)
    z <- stats::rnorm(n)

    y <-
      1 +
      1.7 *
        x +
      0.4 *
        z +
      stats::rnorm(n)

    affected <- 1:6

    x[affected] <-
      x[affected] +
      8

    y[affected] <-
      y[affected] -
      24

    data_case <- data.frame(
      y = y,
      x = x,
      z = z
    )

    set.seed(66002)

    first <- sap_multiscale_test_v2(
      y ~ x + z,
      data = data_case,
      target_pos = 2L,
      k_grid = c(
        1L,
        2L,
        5L,
        6L
      ),
      B_perm = 39L
    )

    set.seed(66002)

    second <- sap_multiscale_test_v2(
      y ~ x + z,
      data = data_case,
      target_pos = 2L,
      k_grid = c(
        1L,
        2L,
        5L,
        6L
      ),
      B_perm = 39L
    )

    testthat::expect_identical(
      first,
      second
    )
  }
)


testthat::test_that(
  "a noiseless exact model does not create a deletion",
  {
    x <- seq(
      -3,
      3,
      length.out = 100L
    )

    data_case <- data.frame(
      x = x,
      y = 2 + 1.5 * x
    )

    set.seed(66003)

    result <- sap_multiscale_test_v2(
      y ~ x,
      data = data_case,
      target_pos = 2L,
      k_grid = c(
        1L,
        2L,
        5L
      ),
      B_perm = 39L
    )

    testthat::expect_false(
      result$detected
    )

    testthat::expect_identical(
      result$selected_k,
      0L
    )

    testthat::expect_length(
      result$robust_cleaning$selected_set,
      0L
    )
  }
)

