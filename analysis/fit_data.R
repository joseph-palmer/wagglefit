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

#' Run my analysis
run_analysis <- function(n = 5, xtol = 0) {
  data <- get_data()
  best_result <- run_mutliple(
    data$foraging_distance,
    n = n,
    verbose_r = TRUE,
    xtol = xtol
  )
  return(best_result)
}
