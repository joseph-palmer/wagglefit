test_scout_dist <- function(x, p, ls, qn, a) {
  truth <- data.frame(
    A = c(
      0.77021, 0.45037,
      0.09216, 0.37465,
      0.25678, 0.05814,
      0., 0.,
      0., 0.94935
    )
  )
  ans <- scout_dist(x, p, ls, qn, a)
  test_that("scout_dist gives expected result", {
    expect_identical(round(truth, 5), round(ans, 5))
  })
}

test_recruit_dist <- function(x, p, ln, qn, a) {
  truth <- data.frame(
    A = c(
      0.98663, 0.69391,
      0., 0.44143,
      0.11115, 0.,
      0., 0.,
      0., 0.51367
    )
  )
  ans <- recruit_dist(x, p, ln, qn, a)
  test_that("recruit_test gives expected results", {
    expect_identical(round(truth, 5), round(ans, 5))
  })
}

test_model_all <- function(x, ...) {
  argvals <- list(...)
  truth <- data.frame(
    A = c(
      1.75684, 1.14428,
      0.09216, 0.81608,
      0.36793, 0.05814,
      0., 0.,
      0., 1.46302
    )
  )
  ans <- model_all(x, ...)
  test_that("model_all gives expected results", {
    expect_identical(round(truth, 5), round(ans, 5))
  })
}

test_model_scout <- function(x, ...) {
  argvals <- list(...)
  truth <- data.frame(
    A = c(
      1.54042, 0.90073,
      0.18432, 0.7493,
      0.51356, 0.11628,
      0., 0.,
      0., 1.8987
    )
  )
  ans <- model_scout(x, ...)
  test_that("model_scout gives expected results", {
    expect_identical(round(truth, 5), round(ans, 5))
  })
}

test_model_recruit <- function(x, ...) {
  argvals <- list(...)
  truth <- data.frame(
    A = c(
      1.97325, 1.38783,
      0., 0.88286,
      0.22231, 0.,
      0., 0.,
      0., 1.02735
    )
  )
  ans <- model_recruit(x, ...)
  test_that("model_recruit gives expected results", {
    expect_identical(round(truth, 5), round(ans, 5))
  })
}

test_model <- function(x, p, ls, ln, qn, a) {
  truth_list <- list(
    "0" = data.frame(
      A = c(
        1.75684, 1.14428,
        0.09216, 0.81608,
        0.36793, 0.05814,
        0., 0.,
        0., 1.46302
      )
    ),
    "1" = data.frame(
      A = c(
        1.54042, 0.90073,
        0.18432, 0.7493,
        0.51356, 0.11628,
        0., 0.,
        0., 1.8987
      )
    ),
    "2" = data.frame(
      A = c(
        1.97325, 1.38783,
        0., 0.88286,
        0.22231, 0.,
        0., 0.,
        0., 1.02735
      )
    )
  )
  ans0 <- model(x, "0", p, ls, ln, qn, a)
  test_that("model with mtype=0 gives expected results", {
    expect_identical(round(truth_list[["0"]], 5), round(ans0, 5))
  })
  ans1 <- model(x, "1", ls, qn, a)
  test_that("model with mtype=1 gives expected results", {
    expect_identical(round(truth_list[["1"]], 5), round(ans1, 5))
  })
  ans2 <- model(x, "2", ln, qn, a)
  test_that("model with mtype=2 gives expected results", {
    expect_identical(round(truth_list[["2"]], 5), round(ans2, 5))
  })
}

test_psudo_model_loglike <- function(x) {
  truth <- -0.2879841981
  test_that("logic of model_loglike works (e.g. sum(log(data.frame)))", {
    expect_identical(truth, round(sum(log(x)), 10))
  })
}

model_tests <- function() {
  x <- data.frame(
    A = c(
      0.2, 0.5,
      1.3, 0.6,
      0.8, 1.5,
      3.8, 8.5,
      3.1, 0.08
    )
  )
  p <- 0.5
  ls <- 1.3
  ln <- 1.3
  qn <- 2.2
  a <- 0.5

  # run through tests
  test_scout_dist(x, p, ls, qn, a)
  test_recruit_dist(x, p, ln, qn, a)
  test_model_all(x, p, ls, ln, qn, a)
  test_model_scout(x, ls, qn, a)
  test_model_recruit(x, ln, qn, a)
  test_model(x, p, ls, ln, qn, a)
  test_psudo_model_loglike(x)
}


# run tests
model_tests()
