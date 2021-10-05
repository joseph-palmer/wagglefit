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
run_analysis <- function(data, model = "all", n = 5, upper = 5, bounds = NULL,
                         xtol = 0, verbose_r = FALSE, verbose = FALSE) {
  best_result <- fit_mutliple(
    data$foraging_distance,
    model = model,
    n = n,
    upper = upper,
    bounds = bounds,
    verbose_r = verbose_r,
    verbose = verbose,
    xtol = xtol
  )
  best_result$data_name <- model
  return(best_result)
}

#' Run on specific site k times and output likelihood space and plots
fit_collective_model_to_data <- function(data, bounds, k = 10, n = 10,
                                         verbose_r = FALSE, verbose = FALSE) {
  # run optimisation k times and return best one
  result_data <- purrr::map_df(
    seq_len(k),
    ~ {
      print(paste("Itteration", .x))
      result <- run_analysis(
        data,
        model = "collective", n = n, upper = 10, bounds = bounds,
        verbose_r = verbose_r, verbose = verbose
      )
      return(
        data.frame(
          ll = result$fmax,
          p = result$est[1],
          bs = result$est[2],
          br = result$est[3],
          as = result$est[4],
          ar = result$est[5]
        )
      )
    }
  )

  # get the best model params
  best_mod <- result_data %>% filter(ll == max(ll))
  if (length(best_mod[, 1] > 1)) {
    best_mod <- best_mod[1, ]
  }
  params <- best_mod %>% select(p, bs, br, as, ar)
  solution <- list(
    fmax = best_mod$ll,
    data_name = "collective",
    est = params %>% as.double()
  )

  # show likelihood space (useful when working out parameter bounds)
  llspace <- map_likelihood_space(
    data$foraging_distance, params,
    model_name = "collective", n = 1000,
    upper = use_upper, bounds = bounds
  )

  return(list(solution = solution, llspace = llspace))
}

#' Run on specific site k times and output likelihood space and plots
fit_individual_model_to_data <- function(data, bounds, k = 10, n = 10,
                                         verbose_r = FALSE, verbose = FALSE) {
  # run optimisation k times and return best one
  result_data <- purrr::map_df(
    seq_len(k),
    ~ {
      print(paste("Itteration", .x))
      result <- run_analysis(
        data,
        model = "individual", n = n, upper = 10, bounds = bounds,
        verbose_r = verbose_r, verbose = verbose
      )
      return(
        data.frame(
          ll = result$fmax,
          bs = result$est[1],
          as = result$est[2]
        )
      )
    }
  )

  # get the best model params
  best_mod <- result_data %>% filter(ll == max(ll))
  if (length(best_mod[, 1] > 1)) {
    best_mod <- best_mod[1, ]
  }
  params <- best_mod %>% select(bs, as)
  solution <- list(
    fmax = best_mod$ll,
    data_name = "individual",
    est = params %>% as.double()
  )

  # show likelihood space (useful when working out parameter bounds)
  llspace <- map_likelihood_space(
    data$foraging_distance, params,
    model_name = "individual", n = 1000,
    upper = use_upper, bounds = bounds
  )

  return(list(solution = solution, llspace = llspace))
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
