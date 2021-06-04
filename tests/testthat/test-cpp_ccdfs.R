test_heaviside <- function() {
  tests <- c(-10, -1, -0.4, -0.001, 0, 0.6, 1, 100)
  expected <- c(0, 0, 0, 0, 1, 1, 1, 1)
  actual <- purrr::map_dbl(
    tests,
    heaviside
  )
  test_that("heaviside works as expected", {
    expect_identical(expected, actual)
  })
}

test_scout_ccdf <- function(x, p, ls, qn, a) {
  expected <- c(
    0.41257543, 0.2264957,
    0.33921855, 0.18375567,
    0.11885913, 0.2264957,
    0.09461152, 0.18375567,
    1.5, 0.23602505
  )
  actual <- purrr::map_dbl(
    x,
    ~ {
      scout_ccdf(.x, p, ls, qn, a, min(x))
    }
  )
  test_that("scout ccdf returns expected results", {
    expect_identical(round(expected, 5), round(actual, 5))
  })
}

test_recruit_ccdf <- function(x, p, ln, qn, a) {
  expected <- c(
    0.41534078, 0.11970665,
    0.30781376, 0.06187754,
    0.0100956, 0.11970665,
    0.00273772, 0.06187754,
    0.5, 0.13444239
  )
  actual <- purrr::map_dbl(
    x,
    ~ {
      recruit_ccdf(.x, p, ln, qn, a, min(x))
    }
  )
  test_that("scout ccdf returns expected results", {
    expect_identical(round(expected, 5), round(actual, 5))
  })
}

test_ccdf_model_all <- function(x, p, ls, ln, qn, a) {
  expected <- c(
    0.82791621, 0.34620235,
    0.64703231, 0.24563321,
    0.12895473, 0.34620235,
    0.09734924, 0.24563321,
    2.00000000, 0.37046744
  )
  result <- rep(0.0, length(x))
  ccdf_model_all(x, result, p, ls, ln, qn, a)
  test_that("ccdf_model_all returns the expected result", {
    expect_identical(round(expected, 5), round(result, 5))
  })
}

test_ccdf_model_scout <- function(x, ls, qn, a) {
  expected <- c(
    0.82515086, 0.45299139,
    0.67843709, 0.36751134,
    0.23771827, 0.45299139,
    0.18922304, 0.36751134,
    2., 0.4720501
  )
  result <- rep(0.0, length(x))
  ccdf_model_scout(x, result, ls, qn, a)
  test_that("ccdf_model_scout returns the expected result", {
    expect_identical(round(expected, 5), round(result, 5))
  })
}

test_ccdf_model_recruit <- function(x, ln, qn, a) {
  expected <- c(
    0.83068156, 0.23941329,
    0.61562751, 0.12375507,
    0.0201912, 0.23941329,
    0.00547545, 0.12375507,
    1., 0.26888477
  )
  result <- rep(0.0, length(x))
  ccdf_model_recruit(x, result, ln, qn, a)
  test_that("ccdf_model_recruit returns the expected result", {
    expect_identical(round(expected, 5), round(result, 5))
  })
}

ccdf_tests <- function() {
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
  test_heaviside()
  test_scout_ccdf(x, p, ls, qn, a)
  test_recruit_ccdf(x, p, ln, qn, a)
  test_ccdf_model_all(x, p, ls, ln, qn, a)
  test_ccdf_model_scout(x, ls, qn, a)
  test_ccdf_model_recruit(x, ln, qn, a)
}

# run tests
ccdf_tests()
