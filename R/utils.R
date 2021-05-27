#' Wrapper for message
#'
#' @description A wrapper for `message` that only prints output when
#' `verbose = TRUE`.
#' @param verbose Logical, defaults to `TRUE`. Should verbose processing
#' messages and warnings be returned.
#' @param ... Additional arguments passed to `message`.
#' @concept utility
message_verbose <- function(verbose = TRUE, ...) {
  if (verbose) {
    message(...)
  }
  return(invisible(NULL))
}

#' Caluclates distance from waggle dance
#'
#' Calculate distance from waggle run duration
#' @description Calculates distance in Km from dance duration, using a linear
#' relationship. Taken from Schürch et al (2019)
#' @param duration double The duration of the waggle run in seconds.
#' @return The distance indicated by the dance
#' @export
#' @examples
#' \dontrun{
#' calc_dist(1.8)
#' }
#'
calc_dist <- function(duration) {
  intercept <- 0.17
  slope <- 1.38
  return(intercept + (slope * duration))
}

#' Draws samples from a normal distribtion trucated between upr and lwr bounds.
#'
#' @description Samples from a normal distribution and then checks it falls
#' within given bounds.
#' @param n integer Number of samples to draw
#' @param mean double Mean for the normal distribution to sample
#' @param sd double Standard deviation for the normal distribution to sample
#' @param lwr double Lower bound samples must be >= to
#' @param upr double Upper bound samples must be <= to
#' @return x double The sampled value
#' @importFrom stats rnorm
#' @export
#' @examples
#' \dontrun{
#' trunc_normal(1, 0, 0.5, -1, 1)
#' }
#'
trunc_normal <- function(n, mean, sd, lwr, upr) {
  continue <- TRUE
  while (continue) {
    x <- rnorm(n = n, mean = mean, sd = sd)
    if (x >= lwr & x <= upr) {
      continue <- FALSE
    }
  }
  return(x)
}
