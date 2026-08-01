# Convert x to a standard score

Convert x to a standard score

## Usage

``` r
x2standard(
  x,
  mu_x = mean(x, na.rm = T),
  sigma_x = stats::sd(x, na.rm = T),
  mu_new = 100,
  sigma_new = 15,
  digits = ifelse(sigma_new == 1, 2, 0),
  rxx = 0.95,
  ci_level = 0.95
)
```

## Arguments

- x:

  a numeric vector

- mu_x:

  mean of current scores

- sigma_x:

  standard deviation of current scores

- mu_new:

  mean of new scores

- sigma_new:

  standard deviation of new scores

- digits:

  rounding digits

## Value

numeric vector

## Examples

``` r
x2standard(13, mu_x = 10, sigma_x = 3)
#> <psycheval::standard_score> num 115
#> Error in loadNamespace(x): there is no package called ‘WJSmisc’
```
