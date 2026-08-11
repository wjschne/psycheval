# Convert x to a standard score

Convert x to a standard score

## Usage

``` r
x2standard(
  x,
  mu_x = mean(x, na.rm = TRUE),
  sigma_x = stats::sd(x, na.rm = TRUE),
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

- rxx:

  reliability coefficient

- ci_level:

  confidence interval percentage

## Value

numeric vector

## Examples

``` r
x2standard(13, mu_x = 10, sigma_x = 3)
#> <psycheval::standard_score> num 115
#>  @ mu                           : num 100
#>  @ sigma                        : num 15
#>  @ rxx                          : num 0.95
#>  @ ci_level                     : num 0.95
#>  @ estimated_true_score         : num 114
#>  @ estimated_true_score_rescaled: num 115
#>  @ standard_error_of_measurement: num 3.35
#>  @ standard_error_of_estimation : num 3.27
#>  @ z                            : num 1.96
#>  @ margin_of_error              : num 6.41
#>  @ ci_lower_bound               : num 108
#>  @ accuracy                     : num 1
#>  @ ci_upper_bound               : num 121
#>  @ ci                           : chr "108—121"
#>  @ percentile                   : chr ".84"
```
