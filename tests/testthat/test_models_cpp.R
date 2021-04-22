test_scout_dist <- function(x, m, p, ls, qn, a) {
  truth <- 0.8000732527680952
  ans <- scout_dist(x, m, p, ls, qn, a)
  test_that("scout_dist gives expected result", {
    expect_equal(truth, ans)
  })
}

test_recruit_dist <- function(x, m, p, ln, qn, a) {
  truth <- 1.0095546946431222
  ans <- recruit_dist(x, m, p, ln, qn, a)
  test_that("recruit_dist gives expected results", {
    expect_identical(round(truth, 5), round(ans, 5))
  })
}

test_loglike_model_all <- function(x, p, ls, ln, qn, a) {
  truth <- -0.48062578001916395
  ans <- loglike_model_all(x, p, ls, ln, qn, a)
  test_that("loglike_model_all gives expected results", {
    expect_equal(truth, ans)
  })
}

test_loglike_model_scout <- function(x, ls, qn, a) {
  truth <- -0.7072428399318523
  ans <- loglike_model_scout(x, ls, qn, a)
  test_that("loglike_model_all gives expected results", {
    expect_equal(truth, ans)
  })
}

test_loglike_model_recruit <- function(x, ln, qn, a) {
  truth <- -1.3118270596288413
  ans <- loglike_model_recruit(x, ln, qn, a)
  test_that("loglike_model_all gives expected results", {
    expect_equal(truth, ans)
  })
}


model_tests <- function() {
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
  test_scout_dist(x[1], min(x), p, ls, qn, a)
  test_recruit_dist(x[1], min(x), p, ln, qn, a)
  test_loglike_model_all(x, p, ls, ln, qn, a)
  test_loglike_model_scout(x, ls, qn, a)
  test_loglike_model_recruit(x, ln, qn, a)
}


# run tests
model_tests()
