% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset_selection.R
\name{subset_selection}
\alias{subset_selection}
\title{subset selection method}
\usage{
subset_selection(pop, var_type, var, method, lambda)
}
\arguments{
\item{pop}{A list. Each list item contains the sample vector of each population.}

\item{var_type}{Factor, can be either "known", for the known variance case, if so specify the variance at parameter "var". "eq_var", for the unknown but equal variance case, and "uneq_var", for the unknown and unequal variance case.}

\item{var}{Numeric, value of variance when variance is known. If unknown variance, ignore this parameter.}

\item{method}{Factor, can be either "base" for the basic subset selection, AS for the Adaptive method, or RAS for the randomised Adaptive method.}

\item{lambda}{Factor, can be either "MSE" which uses the automatic lambda tuning method, or a predetermined number between 0 and 1.}
}
\value{
A numeric vector of R-values, one for each population in the order of the input data.
}
\description{
This function performs the subset selection method, given normally distributed data, to produce R-values for all populations
}
