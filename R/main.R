#' Convert x to a z-score
#'
#' @param x a numeric vector
#' @param mu mean
#' @param sigma standard deviation
#' @export
#'
x2z <- function(x, mu = mean(x, na.rm = T), sigma = stats::sd(x, na.rm = T)) {
  (x - mu) / sigma
}
