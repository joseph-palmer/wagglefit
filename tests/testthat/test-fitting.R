# generate_bounds
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

test_generate_bounds_all <- function() {
  upper <- 5
  p_bnds <- c(0, 1.0)
  ls_bnds <- c(1.0e-6, upper)
  ln_bnds <- c(0, upper)
  q_bnds <- c(1, upper)
  a_bnds <- c(0, upper)
  actual_bounds <- rbind(
    p_bnds, ls_bnds,
    ln_bnds, q_bnds,
    a_bnds
  )
  test_that("generate_bounds_all generates the correct bounds", {
    expect_identical(actual_bounds, generate_bounds_all(upper))
  })
}

test_generate_bounds_scout <- function() {
  upper <- 5
  ls_bnds <- c(1.0e-6, upper)
  q_bnds <- c(1, upper)
  a_bnds <- c(0, upper)
  actual_bounds <- rbind(
    ls_bnds, q_bnds, a_bnds
  )
  test_that("generate_bounds_all generates the correct bounds", {
    expect_identical(actual_bounds, generate_bounds_scout(upper))
  })
}

test_generate_starting_est_all <- function(x) {
  bnds <- generate_bounds_all(5)
  result <- generate_starting_ests_all(x, bounds = bnds)
  test_that("starting paramaters generated are within specified bounds", {
    purrr::walk(
      seq_len(length(result)),
      ~ {
        expect_true((result[[.x]] > bnds[.x, 1]) & (result[[.x]] < bnds[.x, 2]))
      }
    )
  })
}

test_generate_starting_est_scout <- function(x) {
  bnds <- generate_bounds_scout(5)
  result <- generate_starting_ests_scout(x, bounds = bnds)
  test_that("starting paramaters generated are within specified bounds", {
    purrr::walk(
      seq_len(length(result)),
      ~ {
        expect_true((result[[.x]] > bnds[.x, 1]) & (result[[.x]] < bnds[.x, 2]))
      }
    )
  })
}

test_fit <- function(x, p, ls, ln, qn, a) {
  actual_all <- list(
    "fmax" = 1.093483,
    "est" = c(4.246401e-01, 1.540067, 1.754589e-06, 1.819189, 4.494388e-01)
  )
  actual_scout <- list(
    "fmax" = 0.903998,
    "est" = c(0.000001, 2.147276, 1.006386)
  )
  mockery::stub(
    fit, "generate_bounds_all",
    function(upper) {
      cbind(c(0, 1e-6, 0, 0, 0), c(1, 5, 5, 5, 5))
    }
  )
  mockery::stub(
    fit, "generate_bounds_scout",
    function(upper) {
      cbind(c(1e-6, 0, 0), c(5, 5, 5))
    }
  )
  mockery::stub(
    fit, "generate_starting_ests_all",
    function(distance, bounds, verbose_r) {
      c(p, ls, ln, qn, a)
    }
  )
  mockery::stub(
    fit, "generate_starting_ests_scout",
    function(distance, bounds, verbose_r) {
      c(ls, qn, a)
    }
  )
  fit_all <- fit(x, model = "all", xtol = 1e-6)
  fit_scout <- fit(x, model = "scout", xtol = 1e-6)
  test_that("fit all model returns expected results", {
    expect_identical(round(actual_all$fmax, 5), round(fit_all$fmax, 5))
    expect_identical(round(actual_all$est, 5), round(fit_all$est, 5))
  })
  test_that("fit scout model returns expected results", {
    expect_identical(round(actual_scout$fmax, 5), round(fit_scout$fmax, 5))
    expect_identical(round(actual_scout$est, 5), round(fit_scout$est, 5))
  })
  test_that("fit gives error when model requested is not scout or all", {
    expect_error(
      fit(x, model = "madeup"),
      paste(
        "Model madeup is not known. Must be either 'all' or 'scout'"
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
  test_generate_bounds_all()
  test_generate_bounds_scout()
  test_generate_starting_est_all(x)
  test_generate_starting_est_scout(x)
  test_fit_multiple_logic()
  # test_fit (x, p, ls, ln, qn, a)
  # fails on check as 'nlopt_create' not provided by package 'nloptr'. Despite
  # the fact it works locally. Keep off unless you wish to run locally (raise
  # issue later for fix)
}

# run tests
fitting_tests()
