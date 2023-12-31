% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lambda_tuning.R
\name{lambda_tuning}
\alias{lambda_tuning}
\title{Automatic tuning of lambda through MSE and bootstrap}
\usage{
lambda_tuning(R, method, k, S = 100)
}
\arguments{
\item{R}{A numeric vector of R-values}

\item{method}{Factor, method used for obtaining the number of effective populations. Can be "AS" or "RAS".}

\item{k}{Numeric, total number of populations}

\item{S}{Numeric, number of bootstrap replicates, defaults to 100.}
}
\value{
the number of effective populations K
}
\description{
This function performs the AS method to estimate the number of effective populations K
}
