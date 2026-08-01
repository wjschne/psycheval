# Convert W scores to a probabilities

Convert W scores to a probabilities

## Usage

``` r
w2p(w = 500, refw = 500)
```

## Arguments

- w:

  person ability in w-score units

- refw:

  item difficulty in w-score units

## Value

numeric vector of probabilities

## Examples

``` r
w2p(w = 520, refw = 500)
#> [1] 0.9
```
