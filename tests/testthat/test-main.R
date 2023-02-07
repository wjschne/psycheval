library(psycheval)
library(testthat)


test_that("x2standard", {
  expect_equal(x2standard(x = 1, mu_x = 0, sigma_x = 1), 115)
})


test_that("conditional_covariance", {
  v_names <- c("x1", "x2", "y1", "y2")
  s <- diag(4)
  dimnames(s) <- list(v_names, v_names)
  x <- c(x1 = 1, x2 = 2)
  s_nocolnames <- `colnames<-`(s, NULL)
  s_norownames <- `rownames<-`(s, NULL)
  x_nonames <- `names<-`(x, NULL)
  x_all <- c(x1 = 1, x2 = 2, y1 = 0, y2 = 3)
  x_wrongvariable <- c(x, d1 = 8)
  s_not_square <- rbind(s,s)
  s_not_symmetric <- s
  s_not_symmetric[1,2] <- .5

  s_singular <- s
  s_singular[1,2] <- 1
  s_singular[2,1] <- 1

  expect_error(conditional_covariance(x, c(1,2)), "sigma must be a matrix.")
  expect_error(conditional_covariance(x_nonames, s), "x must be a named vector.")
  expect_error(conditional_covariance(x, s_nocolnames), "sigma must have column names.")
  expect_error(conditional_covariance(x, s_norownames), "sigma must have row names.")
  expect_error(conditional_covariance(x, s_not_square), "sigma must be a symmetric square matrix.")
  expect_error(conditional_covariance(x, s_not_symmetric), "sigma must be a symmetric square matrix.")

  expect_error(conditional_covariance(x_all, s), "There are no variables in sigma that are not in x.")
  expect_error(conditional_covariance(x_wrongvariable, s), "The following variables in x are not in sigma: d1")
  expect_error(conditional_covariance(x, s, mu = c(1,2)), "mu and the columns of sigma are of different length.")
  expect_error(conditional_covariance(x, s, mu = c(x1 = 1, x3 = 2)), "The names in mu and the column names of sigma are not the same.")
  expect_error(conditional_covariance(x, s_singular), "The covariance matrix of the predictor variables is not positive definite.")
})

test_that("composite score", {
  x <- c(1, 1)
  R <- matrix(c(1,0, 0, 1), nrow = 2)
  composite <- composite_score(x = x,
                  R = R,
                  mu_x = 0,
                 sigma_x = 1,
                 mu_composite = 0,
                 sigma_composite = 1)
  expect_equal(composite, sqrt(2))
  composite2 <- composite_score(x = x,
                               R = R,
                               mu_x = 0,
                               sigma_x = 1,
                               mu_composite = 0,
                               sigma_composite = 1,
                               w = c(sqrt(3) / 2, .5))




})
