#' Sample Data for SemiMMA
#'
#' A dataset containing effect sizes, standard deviations, and sample sizes
#' for demonstrating multivariate meta-analysis with the SemiMMA package.
#'
#' @format A matrix with 16 rows and 11 columns:
#' \describe{
#'   \item{y1, y2, y3, y4, y5}{Effect sizes for five outcomes.}
#'   \item{s1, s2, s3, s4, s5}{Standard deviations for the effect sizes.}
#'   \item{n}{Sample size for each study.}
#' }
#' @source Simulated data for demonstration purposes.
#' @name data_matrix
#' @docType data
#' @keywords datasets
#' @examples
#' \donttest{
#' data(data_matrix)
#' n <- data_matrix[, ncol(data_matrix)]
#' y <- data_matrix[, seq(1, ncol(data_matrix) - 1, 2)]
#' s <- data_matrix[, seq(2, ncol(data_matrix) - 1, 2)]
#' result <- SemiMMA(y, s, n)
#' }
NULL
