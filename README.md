
<!-- README.md is generated from README.Rmd. Please edit that file -->

# psycheval

<!-- badges: start -->
<!-- badges: end -->

The psycheval package is a set of functions that can be useful in
psychological evaluations. It accompanies the [*Individual
Psychometrics*](https://individual-psychometrics.rbind.io/) online
textbook.

## Installation

You can install the development version of psycheval like so:

``` r
remotes::install_github("wjschne/psycheval")
```

# Functions

## Convert a variable to standard scores

Suppose you have a scaled score of 12 (*μ* = 10, *σ* = 3) that you want
to convert to a standard score (*μ* = 100, *σ* = 15).

``` r
library(psycheval)
x2standard(12, mu_x = 10, sigma_x = 3)
#> [1] 110
```

To convert to a z-score (*μ* = 0, *σ* = 1):

``` r
x2standard(12,
           mu_x = 10, sigma_x = 3,
           mu_new = 0, sigma_new = 1)
#> [1] 0.67
```

To convert to T-score (*μ* = 50, *σ* = 10):

``` r
x2standard(12,
           mu_x = 10, sigma_x = 3,
           mu_new = 50, sigma_new = 10)
#> [1] 57
```
