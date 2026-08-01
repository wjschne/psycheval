# General a multivariate confidence interval for a set of scores

General a multivariate confidence interval for a set of scores

## Usage

``` r
multivariate_ci(x, r_xx, mu, sigma, ci = 0.95, v_names = names(x))
```

## Arguments

- x:

  a vector of scores

- r_xx:

  a vector reliability coefficients

- mu:

  a vector means

- sigma:

  a covariance matrix

- ci:

  confidence level

- v_names:

  a vector of names

## Value

A data frame with the following columns:

- `variable` - Variable names

- `x` - Variable scores

- `r_xx` - Reliability coefficients

- `mu_univariate` - Expected true score estimated from the corresponding
  observed score

- `see_univariate` - Standard error of the estimate computed from the
  corresponding reliability coefficient

- `mu_multivariate` - Expected true score estimated from all observed
  scores

- `see_multivariate` - Standard error of the estimate computed from the
  corresponding reliability coefficient

- `upper_univariate` - upper bound of univariate confidence interval

- `lower_univariate` - lower bound of univariate confidence interval

- `upper_multivariate` - upper bound of multivariate confidence interval

- `lower_multivariate` - lower bound of multivariate confidence interval

## Examples

``` r
# Observed Scores
x <- c(
  vci = 130,
  vsi = 130,
  fri = 70,
  wmi = 130,
  psi = 130
)

# Reliability Coefficients
r_xx <- c(
  vci = .92,
  vsi = .92,
  fri = .93,
  wmi = .92,
  psi = .88
)

# Correlation matrix
R <- ("
  index  vci   vsi   fri   wmi   psi
  vci    1.00  0.59  0.59  0.53  0.30
  vsi    0.59  1.00  0.62  0.50  0.36
  fri    0.59  0.62  1.00  0.53  0.31
  wmi    0.53  0.50  0.53  1.00  0.36
  psi    0.30  0.36  0.31  0.36  1.00") |>
  readr::read_tsv() |>
  tibble::column_to_rownames("index") |>
  as.matrix()
#> Warning: The `file` argument of `read_tsv()` should use `I()` for literal data as of
#> readr 2.2.0.
#>   
#>   # Bad (for example):
#>   read_csv("x,y\n1,2")
#>   
#>   # Good:
#>   read_csv(I("x,y\n1,2"))
#> Rows: 5 Columns: 6
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (1): index
#> dbl (5): vci, vsi, fri, wmi, psi
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

# Covariance matrix
sigma <- R * 15^2

# Population means#'
mu <- rep(100, 5)

mci <- multivariate_ci(
  x = x,
  r_xx = r_xx,
  mu = mu,
  sigma = sigma
)

mci
#>     variable   x r_xx mu_univariate see_univariate mu_multivariate
#> vci      vci 130 0.92         127.6       4.069398       126.69804
#> vsi      vsi 130 0.92         127.6       4.069398       126.11601
#> fri      fri  70 0.93          72.1       3.827205        77.64056
#> wmi      wmi 130 0.92         127.6       4.069398       127.29155
#> psi      psi 130 0.88         126.4       4.874423       127.38539
#>     see_multivariate upper_univariate lower_univariate upper_multivariate
#> vci         3.911922        135.57587        119.62413          134.36526
#> vsi         3.896561        135.57587        119.62413          133.75313
#> fri         3.685405         79.60118         64.59882           84.86382
#> wmi         3.953840        135.57587        119.62413          135.04093
#> psi         4.803025        135.95369        116.84631          136.79915
#>     lower_multivariate
#> vci           119.0308
#> vsi           118.4789
#> fri            70.4173
#> wmi           119.5422
#> psi           117.9716

# Conditional covariance of true score estimates
attr(mci, "conditional_covariance")
#>            vci        vsi        fri        wmi        psi
#> vci 15.3031321  0.7940793  0.6610628  0.6064452  0.1059996
#> vsi  0.7940793 15.1831840  0.8677571  0.3369836  0.5206229
#> fri  0.6610628  0.8677571 13.5822068  0.5011386  0.1153453
#> wmi  0.6064452  0.3369836  0.5011386 15.6328484  0.5570710
#> psi  0.1059996  0.5206229  0.1153453  0.5570710 23.0690477
```
