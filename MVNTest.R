#' Performs a test for mutlivariate normality based on Cramer's method..
#'
#' @param X input dataset
#' @param numproj the number of desired replications. Default is 100,000.
#' @return the p-value of the test.
#' @export
#' @examples
#' str_length(letters)
#' str_length(c("i", "like", "programming", NA))


MVNtest <- function(X, numproj = 100000)
{
  # note that the value returned is the q-value of the test
  p <- ncol(X)
  n <- nrow(X)
  x <- matrix(rnorm(numproj * p), nrow = p)
  # generate 100,000, standard
  # p-variate
  # normal random variables.

  y <- matrix(sqrt(apply(x^2, 2, sum)), nrow = p, ncol = numproj, by = T)

  z <- x / y

  tempdat <- as.matrix(X) %*% z  # this gives rise to a p x numproj matrix
  # called tempdat here

  # perform Shapiro-Wilks' test and calculate individual p-values on each of
  # the numproj observation sets.

  pvals <- as.numeric(matrix(unlist(apply(tempdat, 2, shapiro.test)),
                             ncol=4, by = T)[,2])

  return(min(sort(pvals) * numproj / (1:numproj)))
}
