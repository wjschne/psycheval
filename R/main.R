#' Convert x to a standard score
#'
#' @param x a numeric vector
#' @param mu_x mean of current scores
#' @param sigma_x standard deviation of current scores
#' @param mu_new mean of new scores
#' @param sigma_new standard deviation of new scores
#' @param digits rounding digits
#'
#' @return numeric vector
#' @export
#' @examples
#' x2standard(13, mu_x = 10, sigma_x = 3)
#'
x2standard <- function(x,
                       mu_x = mean(x, na.rm = T),
                       sigma_x = stats::sd(x, na.rm = T),
                       mu_new = 100,
                       sigma_new = 15,
                       digits = ifelse(sigma_new == 1, 2, 0)) {
  round(sigma_new * (x - mu_x) / sigma_x + mu_new, digits)
}

#' Computes covariances of composite scores given a covariance matrix and a weight matrix
#'
#' @export
#' @param Sigma Covariance matrix
#' @param w Weight matrix. Must have the same number of rows as R
#' @param correlation If TRUE, return correlations instead of covariances
#' @examples
#' # Create covariance matrix
#' Sigma <- matrix(0.6, nrow = 5, ncol = 5)
#' diag(Sigma) <- 1
#' # Create weight matrix
#' w <- matrix(0, nrow = 5, ncol = 2)
#' w[1:2,1] <- 1
#' w[3:5,2] <- 1
#' w
#' # covariance matrix of weighted sums
#' composite_covariance(Sigma, w)
composite_covariance <- function(Sigma,w, correlation = FALSE) {
  Sigma_composite <- t(w) %*% Sigma %*% w
  if (correlation) stats::cov2cor(Sigma_composite) else Sigma_composite
}

#' General a multivariate confidence interval for a set of scores
#'
#' @param x a vector of scores
#' @param rxx a vector reliability coefficients
#' @param mu  a vector means
#' @param sigma a covariance matrix
#' @param ci confidence level
#' @param v_names a vector of names
#'
#' @return data.frame
#' @export
#'
#' @examples
#' x_wisc <- c(
#'   vci = 130,
#'   vsi = 130,
#'   fri = 70,
#'   wmi = 130,
#'   psi = 130
#' )
#' rxx_wisc <- c(
#'   vci = .92,
#'   vsi = .92,
#'   fri = .93,
#'   wmi = .92,
#'   psi = .88
#'   )
#' R_wisc <- ("
#'   index	vci 	vsi 	fri 	wmi 	psi
#'   vci  	1.00	0.59	0.59	0.53	0.30
#'   vsi  	0.59	1.00	0.62	0.50	0.36
#'   fri  	0.59	0.62	1.00	0.53	0.31
#'   wmi  	0.53	0.50	0.53	1.00	0.36
#'   psi  	0.30	0.36	0.31	0.36	1.00") |>
#'     readr::read_tsv() |>
#'     tibble::column_to_rownames("index") |>
#'     as.matrix()
#'  multivariate_ci(
#'    x = x_wisc,
#'    rxx = rxx_wisc,
#'    mu = rep(100, 5),
#'    sigma = R_wisc * 225
#'  )
multivariate_ci <- function(x, rxx, mu, sigma, ci = .95, v_names = names(x)) {
  v_observed <- paste0(v_names, "_observed")
  v_true <- paste0(v_names, "_true")
  v_all <- c(v_true, v_observed)
  sigma_true <- `diag<-`(sigma , rxx * diag(sigma))
  sigma_all <- `dimnames<-`(rbind(cbind(sigma_true, sigma_true),
                                  cbind(sigma_true, sigma)),
                            list(v_all,
                                 v_all))
  mu_univariate = rxx * (x - mu) + mu

  lower_p <- (1 - ci) / 2
  upper_p <- 1 - lower_p

  mu_conditional <- mu + sigma_true %*% solve(sigma) %*% (x - mu)
  sigma_conditional <-
    sigma_true - sigma_true %*% solve(sigma) %*% t(sigma_true)
  see_univariate <- sqrt(diag(sigma) * (rxx - rxx ^ 2))
  see_multivariate <- sqrt(diag(sigma_conditional))

  data.frame(
    score = v_names,
    x = x,
    rxx = rxx,
    mu_univariate = mu_univariate,
    see_univariate = see_univariate,
    mu_multivariate = mu_conditional,
    see_multivariate = see_multivariate,
    upper_univariate = stats::qnorm(upper_p, mu_univariate, see_univariate),
    lower_univariate = stats::qnorm(lower_p, mu_univariate, see_univariate),
    upper_multivariate = stats::qnorm(upper_p, mu_conditional, see_multivariate),
    lower_multivariate = stats::qnorm(lower_p, mu_conditional, see_multivariate)
  )


}


