# mRMRe

Parallelized Minimum Redundancy, Maximum Relevance (mRMR) ensemble feature selection.

<!-- badges: start -->
[![R-CMD-check](https://github.com/bhklab/mRMRe/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bhklab/mRMRe/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Package page on CRAN: https://cran.r-project.org/package=mRMRe

## Installation

```r
install.packages("mRMRe")
# development version
remotes::install_github("bhklab/mRMRe")
```

## Usage

```r
library(mRMRe)
data(cgps)
dd <- mRMR.data(data = data.frame(cgps.ic50, cgps.ge))
# classic mRMR: 30 features selected for the first column (IC50)
fit <- mRMR.classic(data = dd, target_indices = 1, feature_count = 30)
solutions(fit)
```
