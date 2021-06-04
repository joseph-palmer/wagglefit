test_inverse_cdf <- function(x) {
  test_that("inverse_cdf works as expected", {
    expected <- tibble::tibble(
      sd = c(0.90, 0.80, 0.60, 0.60, 0.50, 0.50, 0.48, 0.30, 0.20, 0.10),
      prob = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
    )
    expect_identical(round(inverse_ccdf(x), 3), expected)
  })

  test_that("inverse_cdf fails if input is not a double array", {
    x <- tibble::tibble(a = c(0, 4, 2, 6, 8, 3), b = c(0, 5, 1, 8, 3, 9))
    expect_error(inverse_ccdf(x), "x must be a double array")
    expect_error(inverse_ccdf(0.5), "x must be a double array")
  })
}

test_make_ccdf_plot_data <- function(x, param_est, model, npoints) {
  expected <- tibble::tibble(
    x_seq = c(
      0.100, 0.189, 0.278, 0.367, 0.456, 0.544, 0.633, 0.722, 0.811, 0.900
    ),
    cumul_ccdf = c(
      1.000, 0.848, 0.687, 0.534, 0.402, 0.297, 0.219, 0.164, 0.125, 0.097
    )
  )
  actual <- round(make_ccdf_plot_data(x, param_est, model, npoints), 3)
  test_that("make_ccdf_plot_data returns expected results", {
    expect_identical(expected, actual)
  })
}

test_make_base_plot <- function(x) {
  plt <- make_base_plot(x)
  expected_xlab <- "Foraging distance (Km)"
  expected_ylab <- "Ln cumulative probability"
  test_that("make_base_plot returns a ggplot object with the correct labels", {
    expect_true(ggplot2::is.ggplot(plt))
    expect_equal(expected_xlab, plt$labels$x)
    expect_equal(expected_ylab, plt$labels$y)
  })
}

test_make_full_plot <- function(x, model_result_list) {
  plt <- make_full_plot(x, model_result_list)
  print(plot)
}

plotting_tests <- function() {
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
  param_est <- c(p, ls, ln, qn, a)
  example_result <- list(
    est = param_est,
    data_name = "all"
  )

  test_inverse_cdf(x)
  test_make_ccdf_plot_data(x, param_est, "all", 10)
  test_make_base_plot(x)
  test_make_full_plot(x, example_result)
}

# run tests
plotting_tests()
