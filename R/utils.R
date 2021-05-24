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

#' Generates bounds for all parameters
#' @description Generates upper and lower bounds for all parameters
#' @return bounds named list of upper and lower bounds
#' @export
#' @examples
#' \dontrun{
#' generate_bounds()
#' }
#'
generate_bounds <- function() {
  # ...
  bounds <- list(
    "lb" = c(),
    "ub" = c()
  )
  return(bounds)
}

#' Generates starting estimates for all numerical optimisation of all parameters
#' @description Generates starting estimates for all parameters by taking a
#' random value between it's upper and lower bound.
#' @param bounds named list of upper and lower bounds for each parameter
#' @return starting_estimates doubleArray of starting parameter estimates
#' @export
#' @examples
#' \dontrun{
#' generate_starting_estimates(generate_bounds())
#' }
#'
generate_starting_estimates <- function(bounds) {
  return(0)
}

#' fit model to data
#' @description Fits specified model to the data given using MLE.
#' @param distance doubleArray The distance decoded from the waggle dance in
#' meters.
#' @return Not sure yet
#' @export
#' @examples
#' \dontrun{
#' fit(distance)
#' }
#'
fit <- function(distance) {
  bounds <- generate_bounds()
  startest <- generate_starting_estimates(bounds)
  result <- optimise_all(
    distance,
    startest,
    bounds$lower,
    bounds$lower
  )
  return_list <- list(
    "fmax" = result[[1]],
    "est" = result[2:length(result)]
  )
}
