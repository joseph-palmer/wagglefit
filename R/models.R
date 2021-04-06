#' Calculate distribution of scout dances
#' @description Calculates distribution of scout foraging distances
#' @inheritParams model_loglike
#' @param p double proportion of scouts
#' @param ls double scout rate
#' @param qn double quality
#' @param a double alpha value
#' @return The scout distribution
#' @export
#'
scout_dist <- function(x, p, ls, qn, a) {
  m <- min(x)
  maxpart <- (-1 + qn - x * a) / a
  maxpart[maxpart < 0] <- 0
  result <- p * (
    (exp(((-1 + qn + m * a - x * a) * ls) / a) * a * (ls**2) * maxpart) /
      (exp(m * ls) * a - exp(((-1 + qn) * ls) / a) *
        (a + ls - qn * ls + m * a * ls))
  )
  return(result)
}

#' Calculate distribution of recruit dances
#' @description Calculates distribution of recruit foraging distances
#' @inheritParams model_loglike
#' @inheritParams scout_dist
#' @param ln double the recruit rate
#' @importFrom pracma erf erfc
#' @return The scout distribution
#' @export
#'
recruit_dist <- function(x, p, ln, qn, a) {
  m <- min(x)
  maxpart <- -x + ((-1 + qn) / (qn * a))
  maxpart[maxpart < 0] <- 0
  res <- (1 - p) * (
    (2 * exp(-pi * (x**2) * ln) * pi * x * ln * maxpart) /
      (((exp(-(m**2) * pi * ln) * (-1 + qn - m * qn * a)) /
        (qn * a)) + ((-1 + erf(m * sqrt(pi) * sqrt(ln)) +
        erfc((sqrt(pi) * (-1 + qn) * sqrt(ln)) /
          (qn * a))) / (2 * sqrt(ln))))
  )
}

#' Calculate distribution of recruit and scout dances combined
#' @description Calculates superimposed distribution of recruit and scout
#' foraging distances
#' @inheritParams model_loglike
#' @return The scout and recruit model distribution
#' @export
#'
model_all <- function(x, ...) {
  args <- list(...)
  p <- args[[1]]
  ls <- args[[2]]
  ln <- args[[3]]
  qn <- args[[4]]
  a <- args[[5]]
  s_dist <- scout_dist(x, p, ls, qn, a)
  r_dist <- recruit_dist(x, p, ln, qn, a)
  return(s_dist + r_dist)
}

#' Calculate distribution of scout model dances
#' @description Calculates distribution of scout foraging distances for
#' the model
#' @inheritParams model_loglike
#' @return The scout model distribution
#' @export
#'
model_scout <- function(x, ...) {
  args <- list(...)
  ls <- args[[1]]
  qn <- args[[2]]
  a <- args[[3]]
  return(scout_dist(x, 1, ls, qn, a))
}

#' Calculate distribution of recruit model dances
#' @description Calculates distribution of recruit foraging distances for
#' the model
#' @inheritParams model_loglike
#' @return The recruit model distribution
#' @export
#'
model_recruit <- function(x, ...) {
  args <- list(...)
  ln <- args[[1]]
  qn <- args[[2]]
  a <- args[[3]]
  return(recruit_dist(x, 0, ln, qn, a))
}

#' Call desired model
#' @description Calls the requested distribution model (all, scout or recruit)
#' @inheritParams model_loglike
#' @return The requested distribution
#' @export
#'
model <- function(x, mtype, ...) {
  whichmodel <- list(
    "0" = model_all,
    "1" = model_scout,
    "2" = model_recruit
  )
  whichmodel[[mtype]](x, ...)
}

#' Get log-likelihood of requested model
#' @description Calculates log-likelihood of the distribution of the requested
#' models
#' @param x double foraging distances
#' @param mtype char defaults to "0" The model to use
#' 0 = all, 1 = scout, 2 = recruit
#' @param ... parameters for the required model
#' @return The requested distribution log-likelihood
#' @export
#'
model_loglike <- function(x, mtype = "0", ...) {
  return(sum(log(model(x, mtype, ...))))
}
