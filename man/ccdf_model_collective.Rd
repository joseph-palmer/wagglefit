% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ccdf_model_collective}
\alias{ccdf_model_collective}
\title{Model ccdf function for collective model. Stores results in
given array (y)}
\usage{
ccdf_model_collective(x, y, p, bs, br, as, ar)
}
\arguments{
\item{x}{NumericVector foraging distances}

\item{y}{NumericVector storage array for the results}

\item{p}{double Proportion of scouts (0<=p<=1)}

\item{bs}{double Scout rate}

\item{br}{double Recruit rate}

\item{as}{double scout alpha}

\item{ar}{double recruit alpha}
}
\description{
Model ccdf function for collective model. Stores results in
given array (y)
}
