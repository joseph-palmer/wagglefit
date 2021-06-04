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
  expected <- c(0, 1, 2)
  results <- purrr::map_dbl(
    c("all", "scout", "recruit"),
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

test_that("calc_ks works as expected", {
  # add code to test the calc_ks function. Write in once bootstrapped.
  expect_equal(1, 1)
})
