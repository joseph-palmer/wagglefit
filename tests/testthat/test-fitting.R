# get_bounds
test_check_upper_bound <- function() {
  test_that("check_upper_bound gives error if upper < 0", {
    expect_error(check_upper_bound(-1), "upper cannot be < 0")
  })
  test_that("check_upper_bound gives error if upper is not a double", {
    expect_error(check_upper_bound("1"), "upper must be of type 'double'")
  })
  test_that("check_upper_bound gives error if upper is an array", {
    expect_error(check_upper_bound(c(1, 2)), "only pass one value as upper")
  })
}

test_get_bounds_collective <- function() {
  upper <- 5
  p_bnds <- c(0, 1.0)
  bs_bnds <- c(1.0e-6, 10)
  br_bnds <- c(1.0e-6, 10)
  as_bnds <- c(1.0e-12, upper)
  ar_bnds <- c(1.0e-12, upper)
  actual_bounds <- rbind(
    p_bnds, bs_bnds,
    br_bnds, as_bnds,
    ar_bnds
  )
  test_that("get_bounds_collective gets the correct bounds", {
    expect_identical(actual_bounds, get_bounds_collective(upper))
  })
}

test_get_bounds_individual <- function() {
  upper <- 5
  bs_bnds <- c(1.0e-6, upper)
  as_bnds <- c(1.0e-6, upper)
  actual_bounds <- rbind(
    bs_bnds,
    as_bnds
  )
  test_that("get_bounds_collective gets the correct bounds", {
    expect_identical(actual_bounds, get_bounds_individual(upper))
  })
}

test_get_starting_est_collective <- function(x) {
  bnds <- get_bounds_collective(5)
  result <- get_starting_ests_collective(x, bounds = bnds)
  test_that("starting paramaters getd are within specified bounds", {
    purrr::walk(
      seq_len(length(result)),
      ~ {
        expect_true((result[[.x]] > bnds[.x, 1]) & (result[[.x]] < bnds[.x, 2]))
      }
    )
  })
}

test_get_starting_est_individual <- function(x) {
  bnds <- get_bounds_individual(5)
  result <- get_starting_ests_individual(x, bounds = bnds)
  test_that("starting paramaters getd are within specified bounds", {
    purrr::walk(
      seq_len(length(result)),
      ~ {
        expect_true((result[[.x]] > bnds[.x, 1]) & (result[[.x]] < bnds[.x, 2]))
      }
    )
  })
}

test_fit <- function(x, p, ls, ln, qn, a) {
  actual_collective <- list(
    "fmax" = 1.093483,
    "est" = c(4.246401e-01, 1.540067, 1.754589e-06, 1.819189, 4.494388e-01)
  )
  actual_individual <- list(
    "fmax" = 0.903998,
    "est" = c(0.000001, 2.147276, 1.006386)
  )
  mockery::stub(
    fit, "get_bounds_collective",
    function(upper) {
      cbind(c(0, 1e-6, 0, 0, 0), c(1, 5, 5, 5, 5))
    }
  )
  mockery::stub(
    fit, "get_bounds_individual",
    function(upper) {
      cbind(c(1e-6, 0, 0), c(5, 5, 5))
    }
  )
  mockery::stub(
    fit, "get_starting_ests_collective",
    function(distance, bounds, verbose_r) {
      c(p, ls, ln, qn, a)
    }
  )
  mockery::stub(
    fit, "get_starting_ests_individual",
    function(distance, bounds, verbose_r) {
      c(ls, qn, a)
    }
  )
  fit_collective <- fit(x, model = "collective", xtol = 1e-6)
  fit_individual <- fit(x, model = "individual", xtol = 1e-6)
  test_that("fit collective model returns expected results", {
    expect_identical(
      round(actual_collective$fmax, 5), round(fit_collective$fmax, 5)
    )
    expect_identical(
      round(actual_collective$est, 5), round(fit_collective$est, 5)
    )
  })
  test_that("fit individual model returns expected results", {
    expect_identical(
      round(actual_individual$fmax, 5), round(fit_individual$fmax, 5)
    )
    expect_identical(
      round(actual_individual$est, 5), round(fit_individual$est, 5)
    )
  })
  test_that(
    "fit gives error when model requested is not individual or collective", {
    expect_error(
      fit(x, model = "madeup"),
      paste(
        "Model madeup is not known. Must be either 'collective' or 'individual'"
      )
    )
  })
}

test_fit_multiple_logic <- function() {
  result_example <- list(
    list(
      "fmax" = 101,
      "est" = c(0, 0, 0, 0, 0)
    ),
    list(
      "fmax" = 102,
      "est" = c(1, 1, 1, 1, 1)
    ),
    list(
      "fmax" = 1005,
      "est" = c(2, 2, 2, 2, 2)
    )
  )
  results_fmax <- unlist(
    purrr::map(
      result_example, ~ {
        .x[1][[1]]
      }
    )
  )
  max_idx <- which(results_fmax == max(results_fmax))
  test_that("Method used to get the max fit from multiple fits works", {
    expect_identical(result_example[[3]], result_example[[max_idx]])
  })
}

fitting_tests <- function() {
  x <- c(
    0.2, 0.5,
    0.3, 0.6,
    0.8, 0.5,
    0.9, 0.6,
    0.1, 0.48
  )
  p <- 0.5
  ls <- 1.3
  ln <- 1.3
  qn <- 2.2
  a <- 0.5

  # run through model tests
  test_check_upper_bound()
  test_get_bounds_collective()
  test_get_bounds_individual()
  test_get_starting_est_individual(x)
  test_get_starting_est_collective(x)
  test_fit_multiple_logic()
  # test_fit (x, p, ls, ln, qn, a)
  # fails on check as 'nlopt_create' not provided by package 'nloptr'. Despite
  # the fact it works locally. Keep off unless you wish to run locally
}

# run tests
fitting_tests()
