#' Generates bounds for all parameters of the scout recruit superposition model
#'
#' @description Generates upper and lower bounds for all parameters
#' @param upper double The upper bound for the params
#' @return bounds matrix of lower ([,1]) and upper ([,2]) bounds
#' @export
#' @examples
#' \dontrun{
#' generate_bounds_all()
#' }
#'
generate_bounds_all <- function(upper) {
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

#' Generates bounds for all parameters of the scout model
#'
#' @description Generates upper and lower bounds for all parameters
#' @param upper double The upper bound for the params
#' @return bounds matrix of lower ([,1]) and upper ([,2]) bounds
#' @export
#' @examples
#' \dontrun{
#' generate_bounds_all()
#' }
#'
generate_bounds_scout <- function(upper) {
  upper <- as.double(upper)
  ls_bnds <- c(1.0e-6, upper)
  q_bnds <- c(1, upper)
  a_bnds <- c(0, upper)
  bounds <- rbind(
    ls_bnds,
    q_bnds,
    a_bnds
  )
  return(bounds)
}


#' Generates starting estimates for all numerical optimisation of all parameters
#' of the scout recruit superposition model
#'
#' @description Generates starting estimates for all parameters by taking a
#' random value sampled from a normal distribution truncated at the upper and
#' lower bounds. Mean is given as the mid point between the upper and lower
#' bounds and the sd is set as 1.0, except for the p parameter, which has a mean
#' of 0.15 and sd of 0.05 (following estimations that the proportion of scouts
#' in a colony ranges between 10 - 30%).
#' @param distance doubleArray foraging distance to fit.
#' @param bounds matrix of lower ([,1]) and upper ([,2]) bounds
#' @param verbose Bool, to display function evaluations, defaults to FALSE
#' @importFrom purrr map2_dbl
#' @return starting_estimates doubleArray of starting parameter estimates
#' @export
#' @examples
#' \dontrun{
#' generate_starting_estimates_all(distance, generate_bounds())
#' }
#'
generate_starting_ests_all <- function(distance, bounds, verbose = FALSE) {
  message_verbose(verbose, "  -Generating viable starting estimates")
  res <- Inf
  while (is.infinite(res) | is.nan(res)) {
    non_p_bounds <- bounds[2:nrow(bounds), ]
    p_startest <- trunc_normal(
      n = 1, mean = 0.15, sd = 0.05,
      lwr = as.double(bounds[1, 1]), upr = as.double(bounds[1, 2])
    )
    rest_startest <- map2_dbl(
      as.vector(non_p_bounds[, 1]),
      as.vector(non_p_bounds[, 2]),
      ~ {
        trunc_normal(n = 1, mean = (.x + .y) / 2, sd = 1, lwr = .x, upr = .y)
      }
    )
    startest <- c(p_startest, rest_startest)
    res <- loglike_model_all(
      distance,
      startest[1],
      startest[2],
      startest[3],
      startest[4],
      startest[5]
    )
  }
  message_verbose(
    verbose,
    paste(
      c(
        "   -Starting estimates:",
        paste(startest),
        "\n   -Starting function value:",
        res
      ),
      collapse = " "
    )
  )
  return(startest)
}

#' Generates starting estimates for all numerical optimisation of all parameters
#' of the scout model
#'
#' @description Generates starting estimates for all parameters by taking a
#' random value sampled from a normal distribution truncated at the upper and
#' lower bounds. Mean is given as the mid point between the upper and lower
#' bounds and the sd is set as 1.0
#' @param distance doubleArray foraging distance to fit.
#' @param bounds matrix of lower ([,1]) and upper ([,2]) bounds
#' @param verbose Bool, to display function evaluations, defaults to FALSE
#' @importFrom purrr map2_dbl
#' @return starting_estimates doubleArray of starting parameter estimates
#' @export
#' @examples
#' \dontrun{
#' generate_starting_estimates_all(distance, generate_bounds())
#' }
#'
generate_starting_ests_scout <- function(distance, bounds, verbose = FALSE) {
  message_verbose(verbose, "  -Generating viable starting estimates")
  res <- Inf
  while (is.infinite(res) | is.nan(res)) {
    startest <- map2_dbl(
      as.vector(bounds[, 1]),
      as.vector(bounds[, 2]),
      ~ {
        trunc_normal(n = 1, mean = (.x + .y) / 2, sd = 1, lwr = .x, upr = .y)
      }
    )
    res <- loglike_model_scout(
      distance,
      startest[1],
      startest[2],
      startest[3]
    )
  }
  message_verbose(
    verbose,
    paste(
      c(
        "   -Starting estimates:",
        paste(startest),
        "\n   -Starting function value:",
        res
      ),
      collapse = " "
    )
  )
  return(startest)
}

#' fit either the scout and recruit superposition model or the scout model to
#' data
#'
#' @description Fits specified model to the data given using MLE.
#' @param distance doubleArray The distance decoded from the waggle dance.
#' @param model characterArray The model to run (scout or all)
#' @param upper double The upper parameter bound. Defaults to 5.
#' @param iteration integerArray The number of times fit has been called out of
#' the number to be called. Defaults to c(1, 1)
#' @param verbose_r Bool, to display R function evaluations, defaults to FALSE
#' @param verbose_cpp Bool, to display C function evaluations, defaults to FALSE
#' @param xtol double, The absolute tolerance on function value. If 0 (default)
#' then default to nlopt default value.
#' @return Not sure yet
#' @export
#' @examples
#' \dontrun{
#' fit_all(distance)
#' }
#'
fit <- function(distance, model = "all", upper = 5, iteration = c(1, 1),
                verbose_r = FALSE, verbose_cpp = FALSE, xtol = 0) {
  model_function_list <- list(
    "all" = c(generate_bounds_all, generate_starting_ests_all),
    "scout" = c(generate_bounds_scout, generate_starting_ests_scout)
  )
  if (!(tolower(model) %in% names(model_function_list))) {
    stop(
      paste(
        "Model", model, "is not known. Must be either 'all' or 'scout'"
      )
    )
  }
  message_verbose(
    verbose_r,
    paste0(
      "Fit ", model, ": ", iteration[1], "/", iteration[2]
    )
  )
  bounds <- model_function_list[[model]][[1]](upper)
  startest <- model_function_list[[model]][[2]](distance, bounds, verbose_r)
  message_verbose(
    verbose_r,
    "  -Optimising"
  )
  result <- optimise(
    distance,
    startest,
    bounds[, 1],
    bounds[, 2],
    verbose_cpp,
    xtol
  )
  message_verbose(
    verbose_r,
    paste0(
      c(
        "  -Optimisation complete\n",
        "    -fmax =", result[[1]],
        "\n    - est: ",
        paste(
          result[2:length(result)]
        ),
        "\n----"
      ),
      collapse = " "
    )
  )
  return(list(
    "fmax" = result[[1]],
    "est" = result[2:length(result)]
  ))
}

#' Run the full model on data
#'
#' @description Run the full model on given data once.
#' @param upper the upper limmit of the parameter bounds
#' @inheritParams fit
#' @export
#'
run <- function(distance, model = "all", upper = 5, verbose_r = FALSE,
                verbose_cpp = FALSE, xtol = 0, iteration = c(1, 1)) {
  results <- fit(
    distance,
    model = model,
    upper = upper,
    verbose_r = verbose_r,
    verbose_cpp = verbose_cpp,
    xtol = xtol,
    iteration = iteration
  )
  return(results)
}

#' Runs the model fitting on data multiple times
#'
#' @description Run the full model on data multiple times (n) and returns the
#' model with the largest likelihood. This allows for a more comprehensive
#' search of the model sample space.
#' @param n integer The number of models to run
#' @inheritParams fit
#' @export
#'
run_mutliple <- function(distance, model = "all", n = 5, upper = 5,
                         verbose_r = FALSE, verbose_cpp = FALSE, xtol = 0) {
  results <- purrr::map(
    seq_len(n),
    ~ {
      run(
        distance,
        model = model,
        upper = upper, verbose_r = verbose_r,
        verbose_cpp = verbose_cpp, xtol = xtol, iteration = c(.x, n)
      )
    }
  )
  results_fmax <- unlist(
    purrr::map(
      results, ~ {
        .x[1][[1]]
      }
    )
  )
  max_idx <- which(results_fmax == max(results_fmax))
  return(results[max_idx])
}
