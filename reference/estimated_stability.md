# Estimated stability coefficient for IQ

Uses best fitting model of general ability from Table 3 from Breit, M.,
Scherrer, V., Tucker-Drob, E. M., & Preckel, F. (2024). The stability of
cognitive abilities: A meta-analytic review of longitudinal studies.
Psychological Bulletin. https://doi.org/10.1037/bul0000425

## Usage

``` r
estimated_stability(
  age,
  interval,
  different_battery_family = FALSE,
  b0 = 0.716,
  b1 = 0.003,
  b2 = -0.258,
  b3 = 0,
  b4 = -0.095,
  b5 = -0.138,
  b6 = -0.104
)
```

## Arguments

- age:

  age at initial testing

- interval:

  year between tests

- different_battery_family:

  different battery family (e.g., WAIS and Stanford-Binet). Defaults to
  `FALSE`, meaning either the same test or test from the same battery
  family (e.g., WISC and WAIS))

- b0:

  Horizontal asymptote

- b1:

  Age scaling factor

- b2:

  Age growth rate

- b3:

  Interval-age interaction

- b4:

  Interval scaling factor

- b5:

  Interval growth rate

- b6:

  Different battery family effect

## Value

a correlation stability coefficient

## Examples

``` r
estimated_stability(
  age = 12,
  interval = 6,
  different_battery_family = TRUE
)
#> [1] 0.6711221
```
