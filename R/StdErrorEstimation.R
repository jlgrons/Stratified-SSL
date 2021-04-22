# Updated: 2021-04-22

#' Standard error estimation based on influence function expansion.
#'
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param beta Regression parameter estimate to compute the score function.
#' @param residual Estimated residual.
#' @export
#' @return Standard error estimate.
#'

StdErrorEstimation <- function(X_labeled, X_unlabeled, y, beta, residual){

  X_labeled_intercept <- cbind(1, X_labaled)
  X_all_intercept <- cbind(1, rbind(X_labeled, X_unlabeled))

  n_all <- nrow(X_all_intercept)
  n_labeled <- nrow(X_labeled_intercept)

  score <-  t(X_labeled_intercept) %*% (
    X_labeled_intercept * residual^2) / n_labeled

  expit_deriv <- c(ExpitDerivative(X_all_intercept %*% beta))
  # Ask Molei about this.
  expit_deriv[which(is.na(expit_deriv))] <- 0

  information <- t(X_all_intercept) %*% (X_all_intercept * expit_deriv) / n_all
  inverse_information <- solve(information)
  variance_est <-  inverse_information %*% score %*% inverse_information

  return(list(std_error = sqrt(diag(inverse_information) / n_labeled),
              inverse_information = inverse_information))
}



