# Compute composite score

Compute composite score

## Usage

``` r
composite_score(
  x,
  R,
  mu_x = 100,
  sigma_x = 15,
  mu_composite = 100,
  sigma_composite = 15,
  w = NULL
)
```

## Arguments

- x:

  Vector of subtest scores

- R:

  Subtest score correlation matrix

- mu_x:

  Vector of subtest means

- sigma_x:

  Vector of subtest standard deviations

- mu_composite:

  Composite mean

- sigma_composite:

  Composite standard deviation

- w:

  Vector of weights

## Value

composite score

## Examples

``` r
# Subtest scores
x <- c(12, 14)
R <- matrix(c(1, .6, .6, 1), nrow = 2)
composite_score(
  x = x,
  R = R,
  mu_x = 10,
  sigma_x = 3
)
#> [1] 116.7705
```
