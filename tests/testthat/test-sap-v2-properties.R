
repo_root <- Sys.getenv(
  "MIS_SAP_REPO",
  unset = "/root/mis_project_m"
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
  "adaptive grids match the validated representative grids",
  {
    testthat::expect_identical(
      as.integer(
        sap_adaptive_k_grid_v2(
          n = 500L,
          design_columns = 2L
        )
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
      as.integer(
        sap_adaptive_k_grid_v2(
          n = 1000L,
          design_columns = 2L
        )
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
      as.integer(
        sap_adaptive_k_grid_v2(
          n = 2500L,
          design_columns = 2L
        )
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
      as.integer(
        sap_adaptive_k_grid_v2(
          n = 5000L,
          design_columns = 2L
        )
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
  "adaptive grids obey fraction and degrees-of-freedom bounds",
  {
    for (n in 20:1000) {
      grid <- sap_adaptive_k_grid_v2(
        n = n,
        design_columns = 2L
      )

      maximum_allowed <- min(
        floor(
          0.05 *
            n
        ),
        n - 3L
      )

      testthat::expect_true(
        all(
          grid >= 1L
        )
      )

      testthat::expect_true(
        all(
          grid <= maximum_allowed
        )
      )

      testthat::expect_identical(
        as.integer(grid),
        sort(
          unique(
            as.integer(grid)
          )
        )
      )
    }
  }
)


testthat::test_that(
  "small-sample impossibility is reported explicitly",
  {
    testthat::expect_error(
      sap_adaptive_k_grid_v2(
        n = 19L,
        design_columns = 2L
      ),
      "No positive coalition size"
    )
  }
)


testthat::test_that(
  "invalid grid specifications are rejected",
  {
    testthat::expect_error(
      sap_adaptive_k_grid_v2(
        n = 100.5
      ),
      "`n`"
    )

    testthat::expect_error(
      sap_adaptive_k_grid_v2(
        n = 100L,
        design_columns = 0L
      ),
      "`design_columns`"
    )

    testthat::expect_error(
      sap_adaptive_k_grid_v2(
        n = 100L,
        fixed_small = c(
          1,
          2.5
        )
      ),
      "`fixed_small`"
    )

    testthat::expect_error(
      sap_adaptive_k_grid_v2(
        n = 100L,
        proportions = c(
          0,
          0.05
        )
      ),
      "`proportions`"
    )

    testthat::expect_error(
      sap_adaptive_k_grid_v2(
        n = 100L,
        proportions = 0.10,
        max_fraction = 0.05
      ),
      "no larger"
    )
  }
)

