#' Generate data from a multivariate normal
#'
#' @param n Size of data set.
#' @param mu Mean vector.
#' @param Sigma Covariance matrix.
#' @export
#' @return A matrix of multivariate normal data.


CovariateGen <- function(n, mu, Sigma){

  return(mvrnorm(n, mu, Sigma))
}

