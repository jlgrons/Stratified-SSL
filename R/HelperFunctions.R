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



