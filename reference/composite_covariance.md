# Computes covariances of composite scores given a covariance matrix and a weight matrix

Computes covariances of composite scores given a covariance matrix and a
weight matrix

## Usage

``` r
composite_covariance(Sigma, w, correlation = FALSE)
```

## Arguments

- Sigma:

  Covariance matrix

- w:

  Weight matrix. Must have the same number of rows as R

- correlation:

  If TRUE, return correlations instead of covariances

## Examples

``` r
# Create variable names
v_names <- c(paste0("A_", 1:3), paste0("B_", 1:3))
v_composites <- c("A", "B")

# Create covariance matrix
Sigma <- matrix(0.6, nrow = 6, ncol = 6, dimnames = list(v_names, v_names))
diag(Sigma) <- 1

# Create weight matrix
w <- matrix(0, nrow = 6, ncol = 2, dimnames = list(v_names, v_composites))
w[v_names[1:3], "A"] <- 1
w[v_names[4:6], "B"] <- 1
w
#>     A B
#> A_1 1 0
#> A_2 1 0
#> A_3 1 0
#> B_1 0 1
#> B_2 0 1
#> B_3 0 1
# covariance matrix of weighted sums
composite_covariance(Sigma, w)
#>     A   B
#> A 6.6 5.4
#> B 5.4 6.6
```
