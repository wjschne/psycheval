# Given a true score, what is the probability a true score is below a threshold?

Assumes multivaraiate normality of true and observed scores

## Usage

``` r
p_true_score_less_than_threshold(x, threshold, rxx, mu = 100, sigma = 15)
```

## Arguments

- x:

  observed score

- threshold:

  threshold score

- rxx:

  reliability coefficient (must be between 0 and 1, exclusively)

- mu:

  population mean (default = 100)

- sigma:

  population standard deviation (default = 15)

## Value

a probability

## Examples

``` r
# What is the probability that a true score is 70 or less when
# the observed score is 65, the reliability coefficient is .95,
# the population mean is 100, and the population standard
# deviation is 15?
p_true_score_less_than_threshold(
  x = 65,
  threshold = 70,
  rxx = .95,
  mu = 100,
  sigma = 15
)
#> [1] 0.8399214
```
