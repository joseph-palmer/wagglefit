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
  test_that("loglike_model_scout gives expected results", {
    expect_equal(truth, ans)
  })
}

test_loglike_model_recruit <- function(x, ln, qn, a) {
  truth <- -1.3118270596288413
  ans <- loglike_model_recruit(x, ln, qn, a)
  test_that("loglike_model_recruit gives expected results", {
    expect_equal(truth, ans)
  })
}

test_optimise_model <- function(x, p, ls, ln, qn, a) {
  actual_all <- c(
    1.093483, 4.246401e-01, 1.540067, 1.754589e-06, 1.819189, 4.494388e-01
  )
  result_all <- optimise_model(
    x, c(p, ls, ln, qn, a), c(0, 1e-6, 0, 0, 0), c(1, 5, 5, 5, 5),
    xtol = 1e-6, model = 0
  )
  test_that("optimise_model all returns expected results", {
    expect_identical(round(actual_all, 5), round(result_all, 5))
  })
  actual_scout <- c(
    0.903998, 0.000001, 2.147276, 1.006386
  )
  result_scout <- optimise_model(
    x, c(ls, qn, a), c(1e-6, 0, 0), c(5, 5, 5),
    xtol = 1e-6, model = 1
  )
  test_that("optimise_model scout returns expected results", {
    expect_identical(round(actual_scout, 5), round(result_scout, 5))
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
  # test_optimise_model (x, p, ls, ln, qn, a)
  # test_optimise_model fails on check as 'nlopt_create' not provided by
  # package 'nloptr'. Despite the fact it works locally. Keep off unless you
  # wish to run locally (raise issue later for fix)
}

# run tests
model_tests()
