# Updated: 2021-04-19

#' Logit Function
#'
#' Calculates \eqn{log(x) - log(1-x)}.
#'
#' @param x Numeric vector.
#' @return Numeric vector.
Logit <- function(x) {
  return(log(x / (1 - x)))
}

#' Expit Function
#'
#' Calculates \eqn{1 / (1 + exp(-x))}.
#'
#' @param x Numeric vector.
#' @return Numeric vector.
Expit <- function(x) {
  return(1 / (1 + exp(-x)))
}

#' Derivative of the Expit Function
#'
#' Calculates \eqn{exp(x) /  (1 + exp(x))^2}.
#'
#' @param x Numeric vector.
#' @param na_correction Whether to convert NAs to zero.
#' @return Numeric vector.
ExpitDerivative <- function(x, na_correction = F) {
  expit_deriv <- exp(x)/(1+exp(x))^2
  if(na_correction){
    expit_deriv[which(is.na(expit_deriv))] = 0
  }
  return(expit_deriv)
}

#' Truncated Cubic Function
#'
#' Calculates \eqn{((X > x)(X-x))^3}.
#'
#' @param x Numeric vector of interest.
#' @param knot_location Knot location for truncation.
#' @return Numeric vector.
TruncatedCubic <- function(x, knot_location){
  return(((x > knot_location) * (x-knot_location))^3)
}

#' One Hot Encoding
#'
#' @param s Numeric vector of interest containing category for each observation.
#' @return Numeric vectors with the corresponding one hot encoding.
OneHotEncoding <- function(s){

  one_hot_encoding <- c()
  num_categories <- length(unique(s))

  for(k in 1:(num_categories - 1)){
    one_hot_encoding <- cbind(one_hot_encoding, ifelse(s == (k - 1), 1, 0))
  }

  return(one_hot_encoding)
}

#' Basis of Interaction Terms
#'
#' @param X Numeric vector of interest containing category for each observation.
#' @return Numeric vectors with the corresponding one hot encoding.
InteractionBasis <- function(X){

  basis <- X

  p <- length(X[1,])

  for (j in 2:p){
    basis <- cbind(basis, X[,1] * X[,j])
  }

  for (j in 3:p){
    basis <- cbind(basis, X[,2] * X[,j])
  }

  return(basis)
}



