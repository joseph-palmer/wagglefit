suppressMessages(library(dplyr))
library(magrittr)
library(tibble)
suppressMessages(library(lubridate))
library(wagglefit)
suppressMessages(library(rlang))
library(ggplot2)
suppressMessages(library(tidyr))
devtools::load_all()

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
run_analysis <- function(data, model = "all", n = 5, xtol = 0) {
  best_result <- fit_mutliple(
    data$foraging_distance,
    model = model,
    n = n,
    verbose_r = FALSE,
    xtol = xtol
  )
  best_result$data_name <- model
  return(best_result)
}

## run it
run_all_models <- function() {
  n <- 5
  data <- get_data()
  for (i in unique(data$site)) {
    print(paste("Working on", i))
    subdat <- data %>% dplyr::filter(site == i)
    result_recruit <- run_analysis(
      subdat,
      n = n, model = "recruit", xtol = 0
    )
    result_all <- run_analysis(subdat, n = n, model = "all", xtol = 0)
    result_scout <- run_analysis(subdat, n = n, model = "scout", xtol = 0)
    result_list <- list(
      "all" = result_all,
      "scout" = result_scout,
      "recruit" = result_recruit
    )
    saveRDS(object = result_list, file = paste0("analysis/results/", i, ".rds"))
    print(paste(i, "Completed"))
  }
}

## process results
process_results <- function() {
  data <- get_data()
  path <- "analysis/results_version1"
  plots <- purrr::map(
    list.files(path),
    ~ {
      site_name <- gsub(".rds", "", .x)
      subdat <- data %>% dplyr::filter(site == site_name)
      df <- readRDS(paste0(path, "/", .x))
      results_tibble <- make_results_tibble(df) %>%
        rowwise() %>%
        mutate(
          ks_p = calc_ks(
            subdat$foraging_distance, c(p, ls, ln, qn, a),
            data_name, TRUE
          ),
          ks_d = calc_ks(
            subdat$foraging_distance, c(p, ls, ln, qn, a),
            data_name, FALSE
          )
        ) %>%
        ungroup()
      ret <- list(
        site_name = list(
          results_tibble,
          make_full_plot(subdat$foraging_distance, df)
        )
      )
      names(ret) <- site_name
      return(ret)
    }
  )

  lwst_aic_models <- purrr::map_df(
    plots,
    ~ {
      dat <- .x[[1]][[1]] %>%
        slice(which.min(AIC)) %>%
        select(data_name)
      return(list(model_name = dat[[1]]))
    }
  )

  # bar chart of AIC values
  aic_bar <- lwst_aic_models %>%
    group_by(model_name) %>%
    summarise(lowest_AIC = n())

  # combine all the data together
  all_data <- purrr::map_df(
    plots,
    ~ {
      dat <- .x[[1]][[1]]
      dat$data_set <- names(.x)
      return(dat)
    }
  ) %>%
    bind_rows()

  ks_plot_dist <- all_data %>%
    ggplot(aes(x = ks_p)) +
    geom_histogram(fill = "#197bbd") +
    geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
    facet_wrap(~data_name, nrow = 3)

  # count for each model the number where ks says different or not different
  ks_bar <- all_data %>%
    group_by(data_name) %>%
    summarise(
      ndiff = sum(ks_p < 0.05, na.rm = TRUE),
      nsame = sum(ks_p >= 0.05, na.rm = TRUE),
      nerror = sum(is_na(ks_p))
    ) %>%
    rename(model_name = data_name)

  bar_plot <- ks_bar %>%
    full_join(aic_bar, on = model_name) %>%
    select(model_name, ks_not_significant = nsame, lowest_AIC) %>%
    pivot_longer(cols = c(ks_not_significant, lowest_AIC)) %>%
    ggplot(aes(x = model_name, y = value)) +
    geom_bar(aes(fill = name), stat = "identity", position = "dodge") +
    labs(x = "Models", y = "Count")

  multi_plot <- cowplot::plot_grid(
    plots[[6]][[1]][[2]],
    labels = c("A"),
    cowplot::plot_grid(
      bar_plot,
      ks_plot_dist,
      ncol = 2,
      labels = c("B", "C")
    ),
    nrow = 2
  ) %>%
    ggsave("analysis/figures/multi_plot.png", .) # , width = 178, units = "mm")

  # plot all model fits
  just_plots <- purrr::map(
    plots,
    ~ {
      plot <- .x[[1]][[2]] +
        theme(legend.position = "none") +
        labs(x = "", y = "", title = names(.x))
      return(plot)
    }
  )
  plt <- cowplot::plot_grid(plotlist = just_plots, ncol = 4)
  ggsave("analysis/figures/all_fits.png")
}
