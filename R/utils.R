#' Get model number from model name
#'
#' @description Converts the model name to its number for cpp. E.g. 'all' is 0
#' @param model characterArray The model to get the number for
#' @export
#' @concept utility
model_number_from_model <- function(model) {
  model_list <- list(
    "all" = 0,
    "scout" = 1
  )
  if (!(model %in% names(model_list))) {
    stop("Model name not found")
  }
  return(model_list[[model]])
}

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
#' @concept utility
calc_dist <- function(duration) {
  intercept <- 0.17
  slope <- 1.38
  return(intercept + (slope * duration))
}

#' Draws sample from a normal distribtion trucated between upr and lwr bounds.
#'
#' @description Samples from a normal distribution and then checks it falls
#' within given bounds.
#' @param mean double Mean for the normal distribution to sample
#' @param sd double Standard deviation for the normal distribution to sample
#' @param lwr double Lower bound samples must be >= to
#' @param upr double Upper bound samples must be <= to
#' @return x double The sampled value
#' @importFrom stats rnorm
#' @export
#' @examples
#' \dontrun{
#' trunc_normal(0, 0.5, -1, 1)
#' }
#' @concept utility
trunc_normal <- function(mean, sd, lwr, upr) {
  continue <- TRUE
  while (continue) {
    x <- rnorm(n = 1, mean = mean, sd = sd)
    if (x >= lwr & x <= upr) {
      continue <- FALSE
    }
  }
  return(x)
}

#' Calculate AIC values for a model
#'
#' @description Calulcates the AIC values for a given model
#' @param k integer The number of parameters
#' @param ll double The maximum log-likelihood (MLE) of the fitted model
#' @export
#' @concept utility
calc_aic <- function(k, ll) {
  return(2 * k - 2 * ll)
}

#' Run a KS test on the data
#'
#' @description Calulctate the ks test stat and p value for a given model
#' @param x doubleArray Foraging distance data
#' @param param_est doubleArray The parameter estimates for a model
#' @param model_type characterArray the name of the model ran
#' @param pvalue Bool To return the p value (TRUE) or D statistic (FALSE).
#' Defaults to TRUE.
#' @importFrom stats ks.test
#' @importFrom rlang .data
#' @concept utility
calc_ks <- function(x, param_est, model_type, pvalue = TRUE) {
  if (model_type != "all") {
    param_est <- param_est[-1]
  }
  param_est <- param_est[!is.na(param_est)]
  ccdf_data <- inverse_ccdf(x)
  ccdf_model <- model_ccdf(x, param_est, model_number_from_model(model_type))
  if (all(is.nan(ccdf_model))) {
    print(paste(model_type, "has bad results"))
    return(NA)
  }
  ks <- ks.test(ccdf_data$prob, ccdf_model)
  if (pvalue) {
    return(ks$p.value)
  } else {
    return(ks$statistic)
  }
}

#' Run a bootstapped KS test on the data
#'
#' @description Calulctate the ks test stat and p value for a given model
#' @param x doubleArray Foraging distance data
#' @param param_est doubleArray The parameter estimates for a model
#' @param model_type characterArray the name of the model ran
#' @param pvalue Bool To return the p value (TRUE) or D statistic (FALSE).
#' Defaults to TRUE.
#' @importFrom stats ks.test
#' @importFrom rlang .data
#' @importFrom Matching ks.boot
#' @concept utility
calc_ks_boot <- function(x, param_est, model_type, pvalue = TRUE) {
  if (model_type != "all") {
    param_est <- param_est[-1]
  }
  param_est <- param_est[!is.na(param_est)]
  ccdf_data <- inverse_ccdf(x)
  ccdf_model <- model_ccdf(x, param_est, model_number_from_model(model_type))
  if (all(is.nan(ccdf_model))) {
    print(paste(model_type, "has bad results"))
    return(NA)
  }
  ks <- ks.boot(ccdf_data$prob, ccdf_model)
  return(ks$ks.boot.pvalue)
}
