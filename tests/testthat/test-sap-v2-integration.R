
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
  "v2 uses the adaptive grid by default",
  {
    set.seed(55001)

    n <- 500L

    data_case <- data.frame(
      x = stats::rnorm(n)
    )

    data_case$y <-
      1 +
      2 *
        data_case$x +
      stats::rnorm(n)

    set.seed(55002)

    result <- sap_multiscale_test_v2(
      y ~ x,
      data = data_case,
      target_pos = 2L,
      B_perm = 39L,
      use_robust_anchor = FALSE
    )

    testthat::expect_identical(
      result$effective_k_grid,
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
      result$grid_source,
      "bounded_adaptive"
    )
  }
)


testthat::test_that(
  "a non-rejection never creates a cleaning set",
  {
    x <- seq(
      -2,
      2,
      length.out = 100L
    )

    data_case <- data.frame(
      x = x,
      y = 1 + 2 * x
    )

    set.seed(55003)

    result <- sap_multiscale_test_v2(
      y ~ x,
      data = data_case,
      target_pos = 2L,
      k_grid = c(
        1L,
        2L,
        5L
      ),
      B_perm = 39L,
      use_robust_anchor = TRUE
    )

    testthat::expect_false(
      result$detected
    )

    testthat::expect_length(
      result$sensitivity_set,
      0L
    )

    testthat::expect_length(
      result$robust_cleaning$selected_set,
      0L
    )

    testthat::expect_false(
      result$robust_cleaning$
        automatic_deletion_permitted
    )
  }
)

