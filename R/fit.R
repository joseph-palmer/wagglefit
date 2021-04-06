#' Create model starting estimates
#' @description generates random starting values for numical optimisation
#' @inheritParams model_loglike
#'
get_model_starting_estimates <- function(mtype) {
  return(0)
}

#' Create model starting parameter bounds
#' @description generates parameter bounds for numical optimisation
#' @inheritParams model_loglike
#'
get_model_parameter_bounds <- function(mtype) {
  return(0)
}

#' Objective function for wagglefit model
#' @description wrapps `model_loglike` returning the negative log-likelihood
#' @inheritParams model_loglike
#' @return negative log-likelihood for wagglefit model
#'
objective <- function(x, mtype, ...) {
  result <- (-1) * model_loglike(x, mtype, ...)
  ifelse(is.finite(result),
    return(result),
    return(99999)
  )
}

#' Fit wagglefit model
#' @description fits the wagglefit model using numerical optmization
#' @inheritParams model_loglike
#' @importFrom nloptr nloptr
#'
fit_objective <- function(x, mtype) {
  x0 <- get_model_starting_estimates(mtype)
  bounds <- get_model_parameter_bounds(mtype)
  opts <- list(
    "algorithm" = "NLOPT-SLSQP",
    "xtol_rel" = 1.0e-8
  )
  result <- nloptr(
    x0 = x0,
    eval_f = objective,
    opts = opts,
    lb = bounds[["lower"]],
    ub = bounds[["upper"]],
    mtype = mtype,
    ... = ...
  )
}
