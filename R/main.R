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

