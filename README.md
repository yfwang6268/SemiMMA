# SemiMMA: Semiparametric Multivariate Meta-Analysis

**SemiMMA** is an R package for performing Semiparametric multivariate meta-analysis to mitigate outcome reporting bias. It uses a scalable approach with selection models and iterative parameter updates, integrating generalized method of moments (GMM) and inverse propensity weighting (IPW) for robust estimation.

## Dependencies

The `SemiMMA` package requires the following R packages and version:

- **R** (>= 3.5.0)
- **MASS** - For multivariate normal sampling
- **matrixcalc** - For matrix operations
- **Matrix** - For sparse and dense matrix support
- **gmm** - For generalized method of moments estimation
- **SQUAREM** - For squared extrapolation methods

## Installation

You can install `SemiMMA` from GitHub using the `devtools` package:

```R
# Install devtools if you don't have it
install.packages("devtools")

# Install SemiMMA
devtools::install_github("yfwang6268/SemiMMA")
```

## Exmaple

Load the package and analyze multivariate meta-analysis data with effect sizes and standard deviations.

```R
library(SemiMMA)

# Sample data: effect sizes, standard deviations, and sample sizes
# Ref: https://doi.org/10.1002/14651858.CD007044.pub4
data(data_matrix, package = "SemiMMA")

# Prepare data
n <- data_matrix[, ncol(data_matrix)]
y <- data_matrix[, seq(1, ncol(data_matrix) - 1, 2)]
s <- data_matrix[, seq(2, ncol(data_matrix) - 1, 2)]

# Impute missing values
s <- impute_within_study_sd(s, n)

# Run analysis
# This produces a vector of estimated effect sizes and standard errors,
result <- SemiMMA(y, s)
print(result)
```

## Reference

- Setthawong, V., Srisubat, A., Potisat, S., Lojanapiwat, B., & Pattanittum, P. (2023). Extracorporeal shock wave lithotripsy (ESWL) versus percutaneous nephrolithotomy (PCNL) or retrograde intrarenal surgery (RIRS) for kidney stones. *Cochrane Database of Systematic Reviews*, 2023(8), CD007044. [https://doi.org/10.1002/14651858.CD007044.pub4](https://doi.org/10.1002/14651858.CD007044.pub4)



