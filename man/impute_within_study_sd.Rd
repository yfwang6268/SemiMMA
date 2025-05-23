% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SemiMMA.R
\name{impute_within_study_sd}
\alias{impute_within_study_sd}
\title{Impute missing within-study standard deviations}
\usage{
impute_within_study_sd(s, n)
}
\arguments{
\item{s}{A matrix with within-study standard deviations having missing values}

\item{n}{A vector of sample size on each arm}
}
\value{
A matrix with imputed within-study standard deviations no missing values
}
\description{
Impute missing within-study standard deviations
}
