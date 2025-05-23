% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SemiMMA.R
\name{SemiMMA}
\alias{SemiMMA}
\title{Fit Semiparametric Multivariate Meta-Analysis Model}
\usage{
SemiMMA(y, s, n)
}
\arguments{
\item{y}{A matrix with effect sizes}

\item{s}{A matrix with within-study standard deviations}

\item{n}{A vector of sample size on each arm}
}
\value{
A vector contains estimated effect size and the standard deviations
}
\description{
Fits a semiparametric model  using an iterative optimization approach.
}
