test_scout_ccdf <- function(x, m, bs, as) {
  expected <- c(
    0.4204785, 0.0000000,
    0.1269553, 0.0000000,
    0.0000000, 0.0000000,
    0.0000000, 0.0000000,
    1.0000000, 0.0000000
  )
  actual <- purrr::map_dbl(
    x,
    ~ {
      scout_ccdf(.x, m, bs, as)
    }
  )
  test_that("scout ccdf returns expected results", {
    expect_identical(round(expected, 5), round(actual, 5))
  })
}

test_recruit_ccdf <- function(x, m, br, ar) {
  expected <- c(
    0.9502439, 0.6768375,
    0.8748425, 0.5690817,
    0.3669045, 0.6768375,
    0.2807545, 0.5690817,
    1.0000000, 0.6981969
  )
  actual <- purrr::map_dbl(
    x,
    ~ {
      recruit_ccdf(.x, m, br, ar)
    }
  )
  test_that("scout ccdf returns expected results", {
    expect_identical(round(expected, 5), round(actual, 5))
  })
}

test_ccdf_model_collective <- function(x, p, bs, br, as, ar) {
  expected <- c(
    0.6853612, 0.3384187,
    0.5008989, 0.2845408,
    0.1834523, 0.3384187,
    0.1403772, 0.2845408,
    1.0000000, 0.3490985
  )
  result <- rep(0.0, length(x))
  ccdf_model_collective(x, result, p, bs, br, as, ar)
  test_that("ccdf_model_collective returns the expected result", {
    expect_identical(round(expected, 5), round(result, 5))
  })
}

test_ccdf_model_individual <- function(x, bs, as) {
  expected <- c(
    0.4204785, 0.0000000,
    0.1269553, 0.0000000,
    0.0000000, 0.0000000,
    0.0000000, 0.0000000,
    1.0000000, 0.0000000
  )
  result <- rep(0.0, length(x))
  ccdf_model_individual(x, result, bs, as)
  test_that("ccdf_model_individual returns the expected result", {
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
  bs <- 1.3
  br <- 1.3
  as <- 2.2
  ar <- 0.5

  # run through model tests
  test_scout_ccdf(x, min(x), bs, as)
  test_recruit_ccdf(x, min(x), br, ar)
  test_ccdf_model_collective(x, p, bs, br, as, ar)
  test_ccdf_model_individual(x, bs, as)
}

# run tests
ccdf_tests()
