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
x@ci
#> [1] "134—144"
x@percentile
#> Error in loadNamespace(x): there is no package called ‘WJSmisc’
x@estimated_true_score
#> [1] 138.8
```
