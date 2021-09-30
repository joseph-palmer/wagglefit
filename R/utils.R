#' Get model number from model name
#'
#' @description Converts the model name to its number for cpp. E.g. 'all' is 0
#' @param model characterArray The model to get the number for
#' @export
#' @concept utility
model_number_from_model <- function(model) {
  model_list <- list(
    "collective" = 0,
    "individual" = 1
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

#' Run a bootstapped KS test on the data
#'
#' @description Calulctate the ks test stat and p value for a given model
#' @param x doubleArray Foraging distance data
#' @param param_est doubleArray The parameter estimates for a model
#' @param model_type characterArray the name of the model ran
#' @importFrom stats ks.test
#' @importFrom rlang .data
#' @importFrom Matching ks.boot
#' @concept utility
calc_ks_boot <- function(x, param_est, model_type) {
  ccdf_data <- inverse_ccdf(x)
  ccdf_model <- model_ccdf(x, param_est, model_number_from_model(model_type))
  if (all(is.nan(ccdf_model))) {
    print(paste(model_type, "has NA results"))
    return(NA)
  }
  ks <- ks.boot(ccdf_data$prob, ccdf_model)
  return(ks)
}

#' Calculate likelihoods for a specified parameter for the collective model
#'
#' @description Calculate likelihoods for a specified parameter to show the
#' likelihood space for that parameter. This uses the bounds defined for that
#' parameter and gets the specified number of values to show the likelihood
#' space.
#' @param data doubleArray The data to check the parameters against. This should
#' be the data the parameters were fit to.
#' @param var characterArray The name of the variable you want to examine the
#' likelihood space for.
#' @param p double The proportion of scouts.
#' @param bs double The scout rate.
#' @param as double The scout alpha.
#' @param br double The recruit rate.
#' @param ar double The recruit alpha.
#' @param n integer The number of points to sample, defaults to 1000.
#' @param upper integer The upper bound
#' @return tibble The varied parameter and associated likelihood and
#' log-likelihood scores.
#' @importFrom tibble tibble
#' @importFrom purrr map_dbl
#' @concept utility
#'
calc_var_likelihood_collective <- function(
  data, var, params, n = 1000, upper = 5, bounds = NULL
  ) {
  p <- params[[1]]
  bs <- params[[2]]
  br <- params[[3]]
  as <- params[[4]]
  ar <- params[[5]]

  if (is.null(bounds)) {
    bnds <- generate_bounds_all(upper)
  } else {
    bnds <- bounds
  }
  var_bnds <- paste0(var, "_bnds")
  vals <- seq(bnds[var_bnds, 1], bnds[var_bnds, 2], length.out = n)
  result <- map_dbl(
    vals,
    ~ {
      switch(var,
        p = loglike_model_collective(data, .x, bs, br, as, ar),
        bs = loglike_model_collective(data, p, .x, br, as, ar),
        br = loglike_model_collective(data, p, bs, .x, as, ar),
        as = loglike_model_collective(data, p, bs, br, .x, ar),
        ar = loglike_model_collective(data, p, bs, br, as, .x),
        stop(paste("Unknown input", var))
      )
    }
  )
  return(tibble(var = vals, loglike = result, likelihood = exp(result)))
}

#' Calculate likelihoods for a specified parameter for the individual model
#'
#' @description Calculate likelihoods for a specified parameter to show the
#' likelihood space for that parameter. This uses the bounds defined for that
#' parameter and gets the specified number of values to show the likelihood
#' space.
#' @param data doubleArray The data to check the parameters against. This should
#' be the data the parameters were fit to.
#' @param var characterArray The name of the variable you want to examine the
#' likelihood space for.
#' @param p double The proportion of scouts.
#' @param bs double The scout rate.
#' @param as double The scout alpha.
#' @param n integer The number of points to sample, defaults to 1000.
#' @param upper integer The upper bound
#' @return tibble The varied parameter and associated likelihood and
#' log-likelihood scores.
#' @importFrom tibble tibble
#' @importFrom purrr map_dbl
#' @concept utility
#'
calc_var_likelihood_individual <- function(
  data, var, params, n = 1000, upper = 5, bounds = NULL
  ) {
  bs <- params[[1]]
  as <- params[[2]]
  if (is.null(bounds)) {
    bnds <- generate_bounds_all(upper)
  } else {
    bnds <- bounds
  }
  var_bnds <- paste0(var, "_bnds")
  vals <- seq(bnds[var_bnds, 1], bnds[var_bnds, 2], length.out = n)
  result <- map_dbl(
    vals,
    ~ {
      switch(var,
        bs = loglike_model_individual(data, .x, as),
        as = loglike_model_individual(data, bs, .x),
        stop(paste("Unknown input", var))
      )
    }
  )
  return(tibble(var = vals, loglike = result, likelihood = exp(result)))
}

#' Calculate likelihoods for all parameters for a given dataset
#'
#' @description Calculate likelihoods for all parameters to show their
#' likelihood space for a given set of data and optimised model parameters. The
#' idea here is that once you have optimised your parameters you pass them in
#' here with the data and it shows the likelihood smace of each one wrt. the
#' others fixed at their optima.
#' @inheritParams calc_var_likelihood
#' @param params named list or tibble of optimised parameter values
#' @param model_name character array The name of the model to evaluate
#' @return ggplotObj The likelihood plot for each parameter.
#' @importFrom tidyr gather
#' @importFrom purrr map
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap labs
#' @concept utility
#'
map_likelihood_space <- function(data, params, model_name = "collective",
                                n = 1000, upper = 5, bounds = NULL) {
  use_params <- list(
    "individual" = c("bs", "as"),
    "collective" = c("p", "bs", "br", "as", "ar")
  )
  use_model <- list(
    "individual" = calc_var_likelihood_individual,
    "collective" = calc_var_likelihood_collective
  )
  result <- map(
    use_params[[model_name]],
    ~ {
      use_model[[model_name]](
        data, .x, params,
        n = n,
        upper = upper,
        bounds = bounds
      ) %>%
        gather("measurement", "value", -var)
    }
  )
  names(result) <- names(params)
  p_plots <- map(
    names(result),
    ~ {
      result[[.x]] %>%
        filter(value > -1e99) %>%
        ggplot(aes(x = var, y = value)) +
        geom_line() +
        facet_wrap(~measurement, scale = "free_y") +
        labs(x = .x)
    }
  )
  master_plot <- plot_grid(
    plotlist = p_plots,
    ncol = 1,
    labels = names(p_plots)
  )

  return(master_plot)
}
