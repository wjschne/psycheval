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
#' # Create variable names
#' v_names <- c(paste0("A_", 1:3),paste0("B_", 1:3))
#' v_composites <- c("A", "B")
#'
#' # Create covariance matrix
#' Sigma <- matrix(0.6, nrow = 6, ncol = 6, dimnames = list(v_names, v_names))
#' diag(Sigma) <- 1
#'
#' # Create weight matrix
#' w <- matrix(0, nrow = 6, ncol = 2, dimnames = list(v_names, v_composites))
#' w[v_names[1:3],"A"] <- 1
#' w[v_names[4:6],"B"] <- 1
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
#' @param r_xx a vector reliability coefficients
#' @param mu  a vector means
#' @param sigma a covariance matrix
#' @param ci confidence level
#' @param v_names a vector of names
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item `variable` - Variable names
#'   \item `x` - Variable scores
#'   \item `r_xx` - Reliability coefficients
#'   \item `mu_univariate` - Expected true score estimated from the corresponding observed score
#'   \item `see_univariate` - Standard error of the estimate computed from the corresponding reliability coefficient
#'   \item `mu_multivariate` - Expected true score estimated from all observed scores
#'   \item `see_multivariate` - Standard error of the estimate computed from the corresponding reliability coefficient
#'   \item `upper_univariate` - upper bound of univariate confidence interval
#'   \item `lower_univariate` - lower bound of univariate confidence interval
#'   \item `upper_multivariate` - upper bound of multivariate confidence interval
#'   \item `lower_multivariate` - lower bound of multivariate confidence interval
#' }
#'
#' @export
#' @examples
#' # Observed Scores
#' x <- c(
#'   vci = 130,
#'   vsi = 130,
#'   fri = 70,
#'   wmi = 130,
#'   psi = 130
#' )
#'
#' # Reliability Coefficients
#' r_xx <- c(
#'   vci = .92,
#'   vsi = .92,
#'   fri = .93,
#'   wmi = .92,
#'   psi = .88
#'   )
#'
#' # Correlation matrix
#' R <- ("
#'   index	vci 	vsi 	fri 	wmi 	psi
#'   vci  	1.00	0.59	0.59	0.53	0.30
#'   vsi  	0.59	1.00	0.62	0.50	0.36
#'   fri  	0.59	0.62	1.00	0.53	0.31
#'   wmi  	0.53	0.50	0.53	1.00	0.36
#'   psi  	0.30	0.36	0.31	0.36	1.00") |>
#'     readr::read_tsv() |>
#'     tibble::column_to_rownames("index") |>
#'     as.matrix()
#'
#'  # Covariance matrix
#'  sigma <- R * 15 ^ 2
#'
#'  # Population means#'
#'  mu <- rep(100, 5)
#'
#'  mci <- multivariate_ci(
#'    x = x,
#'    r_xx = r_xx,
#'    mu = mu,
#'    sigma = sigma
#'  )
#'
#'  mci
#'
#'  # Conditional covariance of true score estimates
#'  attr(mci, "conditional_covariance")
#'
multivariate_ci <- function(x, r_xx, mu, sigma, ci = .95, v_names = names(x)) {
  v_observed <- paste0(v_names, "_observed")
  v_true <- paste0(v_names, "_true")
  v_all <- c(v_true, v_observed)
  sigma_true <- `diag<-`(sigma , r_xx * diag(sigma))
  sigma_all <- `dimnames<-`(rbind(cbind(sigma_true, sigma_true),
                                  cbind(sigma_true, sigma)),
                            list(v_all,
                                 v_all))
  mu_univariate = r_xx * (x - mu) + mu

  lower_p <- (1 - ci) / 2
  upper_p <- 1 - lower_p

  mu_conditional <- mu + sigma_true %*% solve(sigma) %*% (x - mu)
  sigma_conditional <-
    sigma_true - sigma_true %*% solve(sigma) %*% t(sigma_true)
  see_univariate <- sqrt(diag(sigma) * (r_xx - r_xx ^ 2))
  see_multivariate <- sqrt(diag(sigma_conditional))

  d <- data.frame(
    variable = v_names,
    x = x,
    r_xx = r_xx,
    mu_univariate = mu_univariate,
    see_univariate = see_univariate,
    mu_multivariate = mu_conditional,
    see_multivariate = see_multivariate,
    upper_univariate = stats::qnorm(upper_p, mu_univariate, see_univariate),
    lower_univariate = stats::qnorm(lower_p, mu_univariate, see_univariate),
    upper_multivariate = stats::qnorm(upper_p, mu_conditional, see_multivariate),
    lower_multivariate = stats::qnorm(lower_p, mu_conditional, see_multivariate)
  )

  attr(d, "conditional_covariance") <- sigma_conditional

  return(d)


}



#' Convert W scores to a probabilities
#'
#' @param w person ability in w-score units
#' @param refw item difficulty in w-score units
#'
#' @return numeric vector of probabilities
#' @export
#'
#' @examples
#' w2p(w = 520, refw = 500)
w2p <- function(w = 500, refw = 500) {
  (1 + exp(-(w - refw) / (20 / log(9)))) ^ -1
}


#' Convert logits to W scores
#'
#' @param logit numeric vector of logits
#' @param refw numeric vector of reference W scores
#'
#' @return numeric vector of W scores
#' @export
#'
#' @examples
#' logit2w(2)
#'
logit2w <- function(logit, refw = 500) {
  logit * 20 / log(9) + refw
}



#' Convert W scores to logits
#'
#' @param w numeric vector of W scores
#' @param refw numeric vector of reference W scores
#'
#' @return numeric vector of logits
#' @export
#'
#' @examples
#' w2logit(540)
w2logit <- function(w, refw = 500) {
  log(9) * (w - refw) / 20
}


#' Convert ability (in W scores by default) to relative proficiency index
#'
#' @param x numeric vector of ability scores
#' @param mu numeric vector of ability scores of reference group
#' @param scale number vector of scaling factor. The default value (`log(9) / 20`) assumes that x and mu are W scores.
#' @param criterion numeric proficiency criterion (between 0 and 1, exclusive)
#' @param reverse boolean. If TRUE, the criterion refers to the
#' proficiency of the person instead of the proficiency of the peer
#' group. In other words, the role of the x and mu are reversed.
#' @param interpretation If TRUE, the rpi's print method will provide an interpretation of
#' the relative proficiency.
#'
#' @return numeric
#' @export
#'
#' @examples
#' # What is the probability a person with a W score of 540 can pass
#' # an item that a person with a 500 W score can pass with a
#' # probability of .90?
#' rpi(x = 540, mu = 500, criterion = .9)
#' # Same as above but with an interpretive statement
#' rpi(x = 540, mu = 500, criterion = .9, interpretation = TRUE)
#' # When a person with a W score of 540 has a .9 probability of
#' # passing an item, what is the probability that a person with a W
#' # score of 500 will pass it?
#' rpi(x = 540, mu = 500, criterion = .9, reverse = TRUE, interpretation = TRUE)
rpi <- function(x,
                mu = 500,
                scale = 20 / log(9),
                criterion = .9,
                reverse = FALSE,
                interpretation = FALSE) {
  if (criterion >= 1 | criterion <= 0) stop("criterion must be between 0 and 1, exclusive")
  if (reverse) {
    r <- (1 + exp(log((1 - criterion) / criterion) + (x - mu) / scale )) ^ -1

  } else {
    r <- (1 + exp(-(log(criterion / (1 - criterion)) + (x - mu) / scale))) ^ -1
  }

  class(r) <- c("rpi", class(r))
  attr(r, "criterion") <- criterion
  attr(r, "reverse") <- reverse
  attr(r, "interpretation") <- interpretation
  attr(r, "scale") <- scale
  r

}

#' Format rpi class
#'
#' @param x object with class rpi
#' @param ... additional parameters
#'
#' @return text
#' @noRd
#' @keywords internal
#' @export
format.rpi <- function(x, ...) {

  if (attr(x, "interpretation")) {
    if (attr(x, "reverse")) {
      paste0(
        "When this person has a ",
        WJSmisc::prob_label(attr(x, "criterion"), digits = 2),
        " probability of answering an item correctly, a same-age peer of average ability has a ",
        WJSmisc::prob_label(x, digits = 2),
        " probability of answering it correctly."
      )

    } else {
      paste0(
        "When a same-age peer of average ability has a ",
        WJSmisc::prob_label(attr(x, "criterion"), digits = 2),
        " probability of answering an item correctly, this person has a ",
        WJSmisc::prob_label(x, digits = 2),
        " probability of answering it correctly."
      )

    }


  } else x
}

#' Print object of class rpi
#'
#' @param x object of class rpi
#' @param ...  addtional parameters
#'
#' @return character
#' @noRd
#' @keywords internal
#' @export
print.rpi <- function(x, ...) {
  cat(format(x, ...))
  }










#' Conditional Covariance
#'
#' @param x named numeric vector of predictor scores
#' @param sigma named covariance matrix of predictor and outcome variables
#' @param mu a single numeric mean for all variables or a named vector of means of predictor and outcome variables
#'
#' @return list of conditional means and a covariance matrix
#' \itemize{
#'   \item `mu_conditional` - The means of the outcome variables conditioned on the values of the predictors in vector x.
#'   \item `mu_sigma` - The covariance matrix of the outcome variables conditioned on the values of the predictors in vector x.
#'   \item `descriptives_conditional` - A data frame of means and standard deviations of the outcome variables conditioned on the values of the predictors in vector x.
#'   \item `x` - The predictor scores from the x parameter
#'   \item `sigma` - The unconditional covariance matrix from the sigma parameter
#'   \item `mu` - Anamed vector of unconditional means
#' }
#' @export
#'
#' @examples
#' # Named vector of predictor scores
#' x <- c(A = 1)
#'
#' # Named vector of unconditional means
#' mu <- c(A = 0, B = 0, C = 0)
#'
#' # Unconditional covariance matrix with row and column names
#' sigma <- matrix(c(1, .5, .5,
#'                   .5, 1, .5,
#'                   .5, .5, 1),
#'                 nrow = 3,
#'                 ncol = 3,
#'                 dimnames = list(names(mu),
#'                                 names(mu)))
#'
#' # Conditoinal means and covariance matrix
#' conditional_covariance(x = x, sigma = sigma, mu = mu)
#'
conditional_covariance <- function(x, sigma, mu = 0) {
  if (!("matrix" %in% class(sigma))) stop("sigma must be a matrix.")
  if (is.null(names(x))) stop("x must be a named vector.")
  if (is.null(colnames(sigma))) stop("sigma must have column names.")
  if (is.null(rownames(sigma))) stop("sigma must have row names.")
  if (!all(rownames(sigma) == colnames(sigma))) stop("The row and column names for sigma must be identical.")
  if (!isSymmetric(sigma)) stop("sigma must be a symmetric square matrix.")

  v_names <- colnames(sigma)
  x_names <- names(x)
  y_names <- setdiff(v_names, x_names)
  x_missing <- x_names[is.na(x)]
  x_names <- setdiff(x_names, x_missing)

  if (length(y_names) == 0) stop("There are no variables in sigma that are not in x.")

  if (length(setdiff(x_names, v_names)) > 0) stop(paste0("The following variables in x are not in sigma: ", setdiff(x_names, v_names)))

  if (length(mu) == 1) {
    mu = rep(mu, length(v_names))
    names(mu) <- v_names
  } else if (is.null(names(mu))) {
    if (ncol(sigma) != length(mu)) stop("mu and the columns of sigma are of different length.")
    names(mu) <- v_names
  } else if (!setequal(v_names, names(mu))) {
     stop("The names in mu and the column names of sigma are not the same.")
  } else if (is.null(names(mu))) {
    names(mu) <- x_names
  }


  sigma_x <- sigma[x_names, x_names]
  if (!all(eigen(sigma_x)$values > 0)) stop("The covariance matrix of the predictor variables is not positive definite.")
  invsigma_x <- solve(sigma_x)

  mu_y <- mu[y_names]
  mu_x <- mu[x_names]
  sigma_y <- sigma[y_names, y_names]
  sigma_xy <- sigma[x_names, y_names]
  sigma_yx <- sigma[y_names, x_names]
  mu_y.x <- mu_y + (sigma_yx %*% invsigma_x %*% (x - mu_x))[,1, drop = TRUE]

  sigma_y.x <- sigma_y - sigma_yx %*% invsigma_x %*% sigma_xy

  l <- list(mu_conditional = mu_y.x,
       sigma_conditional = sigma_y.x,
       descriptives_conditional = data.frame(
         construct = names(mu_y.x),
         mu_conditional = mu_y.x,
         sigma_conditional = sqrt(diag(sigma_y.x, names = TRUE))),
       x = x,
       sigma = sigma,
       mu = mu)
  l
}




#' Difference score statistics
#'
#' @param x first score
#' @param y second score
#' @param r_xx reliability of x
#' @param r_yy reliability of y
#' @param r_xy correlation between x and y
#' @param mu population mean of both x and y
#' @param sigma population standard deviation of both x and y
#' @param ci confidence interval of difference score
#' @param tails for significance and prevalance of difference scores
#' @param mu_x population mean of x (defaults to mu)
#' @param mu_y population mean of y (defaults to mu)
#' @param sigma_x population standard deviation of x (defaults to sigma)
#' @param sigma_y population standard deviation of y (defaults to sigma)
#' @return list
#' @export
#'
#' @examples
#' difference_score(
#'   x = 120,
#'   y = 110,
#'   r_xx = .95,
#'   r_yy = .92,
#'   r_xy = .65,
#'   mu = 100,
#'   sigma = 15)
difference_score <- function(x,
                             y,
                             r_xx = .90,
                             r_yy = .90,
                             r_xy = 0,
                             mu = 100,
                             sigma = 15,
                             ci = .95,
                             mu_x = mu,
                             mu_y = mu,
                             sigma_x = sigma,
                             sigma_y = sigma,
                             tails = 2) {

  # Difference score
  d <- x - y
  # variance of difference score
  var_d <- sigma_x ^ 2 + sigma_y ^ 2 - r_xy * sigma_x * sigma_y
  # sd of difference score
  sd_d <- sqrt(var_d)
  # variance of true difference score
  var_d_true <- r_xx * sigma_x ^ 2 + r_yy * sigma_y ^ 2 - r_xy * sigma_x * sigma_y

  # reliability of differenc score
  r_dd <- var_d_true / var_d

  # z-score for confidence interval
  z <- stats::qnorm(1 - (1 - ci) / tails)

  # estimated true difference score
  d_true <- r_dd * (d - (mu_x - mu_y)) + (mu_x - mu_y)
  # lower bound of confidence interval
  d_ci_lb <- d_true - z * sd_d * sqrt(r_dd * (1 - r_dd))
  # upper bound of confidence interval
  d_ci_ub <- d_true + z * sd_d * sqrt(r_dd * (1 - r_dd))

  # standard deviatino of difference scores if true scores are equal
  sd_d_if_true_scores_equal <- sqrt((sigma_x ^ 2) * (1 - r_xx) + (sigma_y ^ 2) * (1 - r_yy))

  # significance-value of difference score
  d_sig <- tails * stats::pnorm(-1 * abs(d) / sd_d_if_true_scores_equal)

  # Prevalence of difference score
  d_prevalence <- tails * stats::pnorm(-1 * abs(d) / sd_d)

  # Return list
  list(
    d = d,
    d_ci_lb = d_ci_lb,
    d_ci_ub = d_ci_ub,
    d_sig = d_sig,
    d_prevalence = d_prevalence,
    sd_d = sd_d,
    r_dd = r_dd,
    x = x,
    y = y,
    r_xx = r_xx,
    r_yy = r_yy,
    r_xy = r_xy,
    mu_x = mu_x,
    mu_y = mu_y,
    sigma_x = sigma_x,
    sigma_y = sigma_y,
    tails = tails,
    ci = ci
  )
}


#' Composite composite score
#'
#' @param x Vector of subtest scores
#' @param R Subtest score correlation matrix
#' @param mu_x Vector of subtest means
#' @param sigma_x Vector of subtest standard deviations
#' @param mu_composite Composite mean
#' @param sigma_composite Composite standard deviation
#'
#' @return composite score
#' @export
#'
#' @examples
#' # Subtest scores
#' x <- c(12, 14)
#' R <- matrix(c(1,.6, .6, 1), nrow = 2)
#' composite_score(x = x,
#'                 R = R,
#'                 mu_x = 10,
#'                sigma_x = 3)
composite_score <- function(x, R,
                            mu_x = 100,
                            sigma_x = 15,
                            mu_composite = 100,
                            sigma_composite = 15) {
  k <- length(x)
  if (length(mu_x) == 1) mu_x <- rep(mu_x, k)
  if (length(mu_x) != length(x)) stop("x and mu_x must be the same length.")
  if (length(sigma_x) == 1) sigma_x <- rep(sigma_x, k)
  if (length(sigma_x) != length(x)) stop("x and mu_x must be the same length.")
  if ((nrow(R) != k) | (ncol(R) != k) | !is.matrix(R)) stop("R must a square matrix with the same size as x.")
  if (length(mu_composite) != 1) stop("mu_composite must be a vector of length 1.")
  if (length(sigma_composite) != 1) stop("sigma_composite must be a vector of length 1.")

  sigma_composite * (sum(x - mu_x) / sqrt(sum(diag(sigma_x) %*% R %*% diag(sigma_x)))) + mu_composite
}



