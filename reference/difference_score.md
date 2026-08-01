# Difference score statistics

Difference score statistics

## Usage

``` r
difference_score(
  x,
  y,
  r_xx = 0.9,
  r_yy = 0.9,
  r_xy = 0,
  mu = 100,
  sigma = 15,
  ci = 0.95,
  mu_x = mu,
  mu_y = mu,
  sigma_x = sigma,
  sigma_y = sigma,
  tails = 2
)
```

## Arguments

- x:

  first score

- y:

  second score

- r_xx:

  reliability of x

- r_yy:

  reliability of y

- r_xy:

  correlation between x and y

- mu:

  population mean of both x and y

- sigma:

  population standard deviation of both x and y

- ci:

  confidence interval of difference score

- mu_x:

  population mean of x (defaults to mu)

- mu_y:

  population mean of y (defaults to mu)

- sigma_x:

  population standard deviation of x (defaults to sigma)

- sigma_y:

  population standard deviation of y (defaults to sigma)

- tails:

  for significance and prevalence of difference scores

## Value

list

## Examples

``` r
difference_score(
  x = 120,
  y = 110,
  r_xx = .95,
  r_yy = .92,
  r_xy = .65,
  mu = 100,
  sigma = 15
)
#> $d
#> [1] 10
#> 
#> $d_ci_lb
#> [1] -1.422461
#> 
#> $d_ci_ub
#> [1] 17.70818
#> 
#> $d_sig
#> [1] 0.06445772
#> 
#> $d_prevalence
#> [1] 0.4255561
#> 
#> $sd_d
#> [1] 12.5499
#> 
#> $r_dd
#> [1] 0.8142857
#> 
#> $x
#> [1] 120
#> 
#> $y
#> [1] 110
#> 
#> $r_xx
#> [1] 0.95
#> 
#> $r_yy
#> [1] 0.92
#> 
#> $r_xy
#> [1] 0.65
#> 
#> $mu_x
#> [1] 100
#> 
#> $mu_y
#> [1] 100
#> 
#> $sigma_x
#> [1] 15
#> 
#> $sigma_y
#> [1] 15
#> 
#> $tails
#> [1] 2
#> 
#> $ci
#> [1] 0.95
#> 
```
