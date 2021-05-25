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
#' @param upper double The upper bound for the params
#' @return bounds matrix of lower ([,1]) and upper ([,2]) bounds
#' @export
#' @examples
#' \dontrun{
#' generate_bounds()
#' }
#'
generate_bounds <- function(upper) {
  upper <- as.double(upper)
  p_bnds <- c(0, 1.0)
  ls_bnds <- c(1.0e-6, upper)
  ln_bnds <- c(0, upper)
  q_bnds <- c(1, upper)
  a_bnds <- c(0, upper)
  bounds <- rbind(
    p_bnds, ls_bnds,
    ln_bnds, q_bnds,
    a_bnds
  )
  return(bounds)
}

#' Generates starting estimates for all numerical optimisation of all parameters
#' @description Generates starting estimates for all parameters by taking a
#' random value between it's upper and lower bound.
#' @param bounds matrix of lower ([,1]) and upper ([,2]) bounds
#' @importFrom purrr map2_dbl
#' @return starting_estimates doubleArray of starting parameter estimates
#' @export
#' @examples
#' \dontrun{
#' generate_starting_estimates(generate_bounds())
#' }
#'
generate_starting_estimates <- function(bounds) {
  startest <- map2_dbl(
    as.vector(bounds[, 1]),
    as.vector(bounds[, 2]),
    ~ {
      runif(n = 1, min = .x, max = .y)
    }
  )
  return(startest)
}

#' fit model to data
#' @description Fits specified model to the data given using MLE.
#' @param distance doubleArray The distance decoded from the waggle dance in
#' meters.
#' @param upper double The upper parameter bound. Defaults to 5.
#' @return Not sure yet
#' @export
#' @examples
#' \dontrun{
#' fit(distance)
#' }
#'
fit <- function(distance, upper = 5) {
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
