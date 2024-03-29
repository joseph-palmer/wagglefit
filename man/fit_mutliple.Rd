% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{fit_mutliple}
\alias{fit_mutliple}
\title{Runs the model fitting (fit) on data multiple times and returns the
paramaters with the highest loglikelihood}
\usage{
fit_mutliple(
  distance,
  model = "collective",
  n = 5,
  upper = 5,
  bounds = NULL,
  verbose_r = FALSE,
  verbose_cpp = FALSE,
  xtol = 0
)
}
\arguments{
\item{distance}{doubleArray The distance decoded from the waggle dance.}

\item{model}{characterArray The model to run (collective or individual).
Defaults to collective.}

\item{n}{integer The number of models to run}

\item{upper}{double The upper parameter bound. Defaults to 5.}

\item{bounds}{named list user supplied bounds. Defaults to NULL. If null the
bounds are automatically generated.}

\item{verbose_r}{Bool, to display R function evaluations, defaults to FALSE}

\item{verbose_cpp}{Bool, to display C function evaluations, defaults to FALSE}

\item{xtol}{double, The absolute tolerance on function value. If 0 (default)
then default to nlopt default value.}
}
\description{
Run the model on data multiple times (n) and returns the
model with the largest likelihood. This allows for a more comprehensive
search of the model sample space.
}
