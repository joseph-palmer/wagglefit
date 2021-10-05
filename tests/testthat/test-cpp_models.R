test_scout_dist <- function(x, m, bs, as) {
  truth <- 4.15163
  ans <- scout_dist(x, m, bs, as)
  test_that("scout_dist gives expected result", {
    expect_equal(truth, round(ans, 5))
  })
}

test_recruit_dist <- function(x, m, br, ar) {
  truth <- 0.63748
  ans <- recruit_dist(x, m, br, ar)
  test_that("recruit_dist gives expected results", {
    expect_identical(round(truth, 5), round(ans, 5))
  })
}

test_loglike_model_collective <- function(x, p, bs, br, as, ar) {
  truth <- -2.233698
  ans <- loglike_model_collective(x, p, bs, br, as, ar)
  test_that("loglike_model_collective gives expected results", {
    expect_equal(truth, round(ans, 6))
  })
}

test_loglike_model_individual <- function(x, br, ar) {
  truth <- -0.83449
  ans <- loglike_model_individual(x, br, ar)
  test_that("loglike_model_individual gives expected results", {
    expect_equal(truth, round(ans, 5))
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
  bs <- 1.3
  br <- 1.3
  as <- 2.2
  ar <- 0.5

  # run through model tests
  test_scout_dist(x[1], min(x), bs, as)
  test_recruit_dist(x[1], min(x), br, ar)
  test_loglike_model_collective(x, p, bs, br, as, ar)
  test_loglike_model_individual(x, br, ar)
  # test_optimise_model (x, p, ls, ln, qn, a)
  # test_optimise_model fails on check as 'nlopt_create' not provided by
  # package 'nloptr'. Despite the fact it works locally. Uncomment to run test
  # locally
}

# run tests
model_tests()
