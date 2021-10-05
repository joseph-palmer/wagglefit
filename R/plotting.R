
#' Calculates the ccdf for given data
#'
#' @description Orders the data smallest to largest and calculates the
#' probability of sampling a value greater than or equal to a given value.
#' @param x The data, e.g. foraging distance
#' @return tibble of the sorted foraging distances and their cumulative
#' probability
#' @importFrom purrr map_dbl
#' @importFrom tibble tibble
#' @export
#'
inverse_ccdf <- function(x) {
  if (!is.double(x) | length(x) < 2) {
    stop("x must be a double array")
  }
  sd <- -sort(-x)
  prob <- map_dbl(
    seq_len(length(x)),
    ~ {
      .x / length(x)
    }
  )
  return(tibble(sd, prob))
}

#' Creates data for plots
#'
#' @description Using the optimised paramater estimates, creates the data for
#' the ccdf plots
#' @param x doubleArray The foraging distance to plot
#' @param param_est doubleArray The paramater estimates of the optimisation
#' @param model characterArray The name of the model used
#' @param npoints integer The number of points to plot the predicted data for.
#' Defaults to 100.
#' @return tibble of model data and cumulative probability
#' @importFrom tibble tibble
#' @export
#'
make_ccdf_plot_data <- function(x, param_est, model, npoints = 100) {
  model_n <- model_number_from_model(model)
  x_seq <- seq(min(x), max(x), length.out = npoints)
  cumul_ccdf <- model_ccdf(x_seq, param_est, model_n)
  cumul_ccdf[[1]] <- 1
  cumul <- tibble(x_seq, cumul_ccdf)
  return(cumul)
}

#' Creates ccdf plot of data
#'
#' @description Creates a ccdf plot for the data provided
#' @param x doubleArray The foraging distance to plot
#' @param logit Bool to log the probabilities or not.
#' @return ggplot plot of cumulative probability for foraging distances in Km
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line theme_set theme
#' theme_classic element_text labs
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
make_base_plot <- function(x, logit = TRUE) {
  theme_set(
    theme_classic() +
      theme(
        text = element_text(family = "URWHelvetica", size = 42)
      )
  )
  if (logit) {
    plt <- inverse_ccdf(x) %>%
      ggplot(aes(x = .data$sd, y = log(.data$prob))) +
      geom_point() +
      labs(x = "Foraging distance (Km)", y = "Ln cumulative probability")
  } else {
    plt <- inverse_ccdf(x) %>%
      ggplot(aes(x = .data$sd, y = .data$prob)) +
      geom_point() +
      labs(x = "Foraging distance (Km)", y = "Ln cumulative probability")
  }
  return(plt)
}

#' Creates plots with model fits from optimised paramaters
#'
#' @description Using the optimised paramater estimates, creates ccdf plots
#' @param x The foraging distance to plot
#' @param model_result_list namedlist The return of `fit`, a list of parameter
#' estimates ($est) and the data name ($data_name).
#' @param logit Bool To log the data or not. Defaults to TRUE
#' @param subplot_coords doubleArray Coordinates for location of inner histogram
#' @return ggplot plot of cumulative probability of foraging distances in Km
#' along with the model fits to this data.
#' @importFrom ggplot2 ggplot aes geom_histogram geom_line theme_set theme
#' theme_classic element_text element_blank element_rect annotation_custom
#' ggplotGrob
#' @importFrom tibble tibble
#' @importFrom purrr map map_df
#' @importFrom rlang .data
#' @export
#'
make_full_plot <- function(x, model_result_list,
                           logit = TRUE,
                           subplot_coords = c(0.5, 2.5, -7.8, -2.5)) {
  cdf_data <- map(
    model_result_list,
    ~ {
      result <- make_ccdf_plot_data(x, .x$est, model = .x$data_name)
      if(any((is.nan(result$cumul_ccdf)))) {
        message("The combination of parameters produces nan values in the ccdf")
      }
      return(result)
    }
  )
  df <- map_df(
    cdf_data,
    I,
    .id = "Model"
  )
  plt <- make_base_plot(x, TRUE) +
    geom_line(
      data = df,
      aes(
        x = .data$x_seq,
        y = log(.data$cumul_ccdf),
        colour = .data$Model
      )
    ) +
    theme(
      legend.position = "none"
    )

  histoplot <- as_tibble(x) %>%
    ggplot(aes(x = .data$value)) +
    geom_histogram(
      bins = 100,
      binwidth = (max(x) - min(x)) / 20,
      col = "white"
    ) +
    theme(
      axis.title.x = element_blank(),
      legend.position = "none",
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", colour = NA)
    )

  full_plot <- plt + annotation_custom(
    ggplotGrob(histoplot),
    xmin = subplot_coords[[1]],
    xmax = subplot_coords[[2]],
    ymin = subplot_coords[[3]],
    ymax = subplot_coords[[4]]
  )

  return(full_plot)
}

#' Converts results list to tibble
#'
#' @description Converts the results list from running the models into a wide
#' format tibble for each model
#' @param result namedList the results list
#' @return tibble of model results
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows as_tibble mutate
#' @importFrom tidyr pivot_wider
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
make_results_tibble <- function(result) {
  . <- NULL
  collective_tibble <- as_tibble(result$collective)
  collective_tibble$parameter <- c("p", "bs", "br", "as", "ar")
  individual_tibble <- as_tibble(result$individual)
  individual_tibble$parameter <- c("bs", "as")
  result_tibble <- bind_rows(
    collective_tibble,
    individual_tibble
  )
  result_tibble <- result_tibble %>%
    pivot_wider(names_from = .data$parameter, values_from = .data$est) %>%
    mutate(
      AIC = calc_aic(rowSums(!(is.na(.))) - 2, .data$fmax)
    ) %>%
    mutate(
      p = ifelse(
        .data$data_name == "individual",
        1,
        .data$p
      )
    )
  return(result_tibble)
}
