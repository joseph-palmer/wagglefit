library(dplyr)
library(magrittr)
library(tibble)
library(lubridate)
library(wagglefit)

#' Get waggle dance data from file and calculate foraging distance
#' @descriptionFor For specific data in the 'data' directory read it in,
#' calculate foraging distance and get just the required columns
#' @importFrom rlang .data
#' @export
#'
get_data <- function() {
  full_data_path <- "analysis/data/FullHBForagingData.csv"
  data <- tibble(read.csv(full_data_path))
  data <- data %>%
    mutate(
      date = ymd(.data$date),
      foraging_distance = calc_dist(.data$duration.seconds)
    ) %>%
    select(
      date,
      site,
      foraging_distance
    )
  return(data)
}

#' Run the full model on some data
#' @description Runs the model fitting on given data
#' @param duration double The duration of the waggle run in seconds.
#' @export
#'
run <- function(data) {
  data <- get_data()
  head(data)
  results <- fit(data$foraging_distance, upper = 5)
  return(results)
}

### run process
run()
