# create standard scores

computes confidence intervals, percentiles

## Usage

``` r
standard_score(
  .data = numeric(0),
  mu = 100L,
  sigma = 15,
  rxx = 0.9,
  ci_level = 0.95,
  accuracy = 1
)
```

## Arguments

- .data:

  obtained score

- mu:

  population mean

- sigma:

  population standard deviation

- rxx:

  reliability coefficient

- ci_level:

  confidence level

- accuracy:

  rounding accuracy

## Slots

- `estimated_true_score`:

  estimated true score

- `standard_error_of_measurement`:

  standard error of measurement

- `standard_error_of_estimation`:

  standard error of estimation

- `margin_of_error`:

  margin of error

- `z`:

  z-score associated with `ci_level`

- `ci_lower_bound`:

  lower bound of CI

- `ci_upper_bound`:

  upper bound of CI

- `ci`:

  CI as a range

- `percentile`:

  percentile

## Examples

``` r
x <- standard_score(140, rxx = .97)
x
#> <psycheval::standard_score> num 140
#>  @ mu                           : int 100
#>  @ sigma                        : num 15
#>  @ rxx                          : num 0.97
#>  @ ci_level                     : num 0.95
#>  @ estimated_true_score         : num 139
#>  @ estimated_true_score_rescaled: num 139
#>  @ standard_error_of_measurement: num 2.6
#>  @ standard_error_of_estimation : num 2.56
#>  @ z                            : num 1.96
#>  @ margin_of_error              : num 5.02
#>  @ ci_lower_bound               : num 134
#>  @ accuracy                     : num 1
#>  @ ci_upper_bound               : num 144
#>  @ ci                           : chr "134—144"
#>  @ percentile                   : chr ".996"
x@ci
#> [1] "134—144"
x@percentile
#> [1] ".996"
x@estimated_true_score
#> [1] 138.8
```
