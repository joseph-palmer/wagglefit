% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{make_base_plot}
\alias{make_base_plot}
\title{Creates ccdf plot of data}
\usage{
make_base_plot(x, logit = TRUE)
}
\arguments{
\item{x}{doubleArray The foraging distance to plot}

\item{logit}{Bool to log the probabilities or not.}
}
\value{
ggplot plot of cumulative probability for foraging distances in Km
}
\description{
Creates a ccdf plot for the data provided
}
