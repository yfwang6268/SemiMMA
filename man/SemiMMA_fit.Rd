% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SemiMMA.R
\name{SemiMMA_fit}
\alias{SemiMMA_fit}
\title{Fit Semiparametric Multivariate Meta-Analysis Model}
\usage{
SemiMMA_fit(y, s, J)
}
\arguments{
\item{y}{A matrix with effect sizes}

\item{s}{A matrix with standard deviations}
}
\value{
A matrix with estimates of effect sizes and standard deviations per row
}
\description{
Fits a semiparametric model  using an iterative optimization approach.
}
