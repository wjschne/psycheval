# Conditional Covariance

Conditional Covariance

## Usage

``` r
conditional_covariance(x, sigma, mu = 0)
```

## Arguments

- x:

  named numeric vector of predictor scores

- sigma:

  named covariance matrix of predictor and outcome variables

- mu:

  a single numeric mean for all variables or a named vector of means of
  predictor and outcome variables

## Value

list of conditional means and a covariance matrix

- `mu_conditional` - The means of the outcome variables conditioned on
  the values of the predictors in vector x.

- `mu_sigma` - The covariance matrix of the outcome variables
  conditioned on the values of the predictors in vector x.

- `descriptives_conditional` - A data frame of means and standard
  deviations of the outcome variables conditioned on the values of the
  predictors in vector x.

- `x` - The predictor scores from the x parameter

- `sigma` - The unconditional covariance matrix from the sigma parameter

- `mu` - Anamed vector of unconditional means

- `beta` - Regression Coefficients for finding mu_conditional

## Examples

``` r
# Named vector of predictor scores
x <- c(A = 1)

# Named vector of unconditional means
mu <- c(A = 0, B = 0, C = 0)

# Unconditional covariance matrix with row and column names
sigma <- matrix(
  c(
    1, .5, .5,
    .5, 1, .5,
    .5, .5, 1
  ),
  nrow = 3,
  ncol = 3,
  dimnames = list(
    names(mu),
    names(mu)
  )
)

# Conditoinal means and covariance matrix
conditional_covariance(x = x, sigma = sigma, mu = mu)
#> $mu_conditional
#>   B   C 
#> 0.5 0.5 
#> 
#> $sigma_conditional
#>      B    C
#> B 0.75 0.25
#> C 0.25 0.75
#> 
#> $descriptives_conditional
#>   construct mu_conditional sigma_conditional
#> B         B            0.5         0.8660254
#> C         C            0.5         0.8660254
#> 
#> $x
#> A 
#> 1 
#> 
#> $sigma
#>     A   B   C
#> A 1.0 0.5 0.5
#> B 0.5 1.0 0.5
#> C 0.5 0.5 1.0
#> 
#> $mu
#> A B C 
#> 0 0 0 
#> 
#> $beta
#>      [,1]
#> [1,]  0.5
#> [2,]  0.5
#> 
```
