test_that("calc_dist gives correct answer with double", {
  expect_equal(calc_dist(1.8), 2.654)
})

test_that("cacl_dist gives correct answer with data frame to 3dp", {
  duration <- data.frame(A = c(0.2, 0.6, 0.8, 1.3, 1.6, 4))
  truth <- data.frame(A = c(0.446, 0.998, 1.274, 1.964, 2.378, 5.690))
  distance <- calc_dist(duration)
  expect_identical(round(truth, 3), round(distance, 3))
})

test_that("trunc_normal gives expected results", {
  check <- purrr::map(
    seq_len(10000),
    ~ {
      trunc_normal(mean = 0.2, sd = 0.5, lwr = 0, upr = 1)
    }
  )
  expect_true(
    all((check < 1 & check > 0))
  )
})

test_that("message_verbose returns a message when asked", {
  expect_equal(
    capture.output(message_verbose(TRUE, "hi there"), type = "message"),
    "hi there"
  )
})

test_that("message_verbose does not return a message when asked", {
  l <- suppressMessages(
    capture.output(message_verbose(FALSE, "hi there"), type = "message")
  )
  expect_equal(l, character(0))
})

test_that("model_number_from_model returns correct models", {
  expected <- c(0, 1)
  results <- purrr::map_dbl(
    c("collective", "individual"),
    model_number_from_model
  )
  expect_identical(expected, results)
})

test_that("model_number_from_model returns error if model name is not known", {
  expect_error(
    model_number_from_model("amadeupmodel"),
    "Model name not found"
  )
})

test_that("calc_aic works as expected", {
  expected <- 314
  result <- calc_aic(2, -155)
  expect_equal(result, expected)
})

expect_collective_dstat <- 0.4
expect_individual_dstat <- 0.7
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
param_est_collective <- c(p, bs, br, as, ar)
param_est_individual <- c(bs, as)
collective_result <- calc_ks_boot(x, param_est_collective, "collective")
individual_result <- calc_ks_boot(x, param_est_individual, "individual")
test_that("calc_ks_boot collective gives correct D stat", {
  expect_equal(
    as.numeric(collective_result$ks$statistic), expect_collective_dstat
  )
})
test_that("calc_ks_boot individual gives correct D stat", {
  expect_equal(
    as.numeric(individual_result$ks$statistic), expect_individual_dstat
  )
})
