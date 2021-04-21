test_scout_dist <- function(x, m, p, ls, qn, a) {
  truth <- 0.9696085
  ans <- scout_dist(x, m, p, ls, qn, a)
  test_that("scout_dist gives expected result", {
    expect_identical(round(truth, 5), round(ans, 5))
  })
}

test_recruit_dist <- function(x, m, p, ln, qn, a) {
  truth <- 1.215333
  ans <- recruit_dist(x, m, p, ln, qn, a)
  test_that("recruit_test gives expected results", {
    expect_identical(round(truth, 5), round(ans, 5))
  })
}



model_tests <- function() {
  x <- data.frame(
    A = c(
      0.2, 0.5,
      1.3, 0.6,
      0.8, 1.5,
      3.8, 8.5,
      3.1, 0.48
    )
  )
  p <- 0.5
  ls <- 1.3
  ln <- 1.3
  qn <- 2.2
  a <- 0.5

  # run through tests
  test_scout_dist(x[1, 1], min(x), p, ls, qn, a)
  test_recruit_dist(x[1, 1], min(x), p, ln, qn, a)
}


# run tests
model_tests()
