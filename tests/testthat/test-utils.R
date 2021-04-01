test_that("calc_dist gives correct answer with double", {
  expect_equal(calc_dist(1.8), 2.654)
})
test_that("cacl_dist gives correct answer with data frame to 3dp", {
  duration <- data.frame(A = c(0.2, 0.6, 0.8, 1.3, 1.6, 4))
  truth <- data.frame(A = c(0.446, 0.998, 1.274, 1.964, 2.378, 5.690))
  distance <- calc_dist(duration)
  expect_identical(round(truth, 3), round(distance, 3))
})
