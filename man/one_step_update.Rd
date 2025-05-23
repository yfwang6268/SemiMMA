% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SemiMMA.R
\name{one_step_update}
\alias{one_step_update}
\title{Perform One-Step Parameter Update}
\usage{
one_step_update(outcome_parameter, y, s)
}
\arguments{
\item{outcome_parameter}{A vector with previous parameters}

\item{y}{A matrix with effect sizes}

\item{s}{A matrix with standard deviations}
}
\value{
A vector with updated parameters
}
\description{
Updates the parameter in a single iteration
}
