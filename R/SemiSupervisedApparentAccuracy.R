# Updated: 2021-04-22

#' Apparent estimates for brier score and overall misclassification rate.
#'
#' @param basis_labeled Basis matrix for labeled data set.
#' @param bais_unlabeled Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of weights.
#' @param resamp_weight Numeric vector of resampling weights.
#' @param threshold Threshold for overall misclassification rate.
#' @export
#' @return resamp_weightector containing regression coefficients.
#'

SemiSupervisedApparentAccuracy <- function(basis_labeled, basis_unlabeled,
                                           X_labeled, X_unlabaled, y, beta_SSL,
                                           beta_imp, samp_prob,
                                           resamp_weight = NULL,
                                           threshold = 0.5){

  if(is.null(resamp_weight)){resamp_weight <- rep(1, length(y))}
  weight <- resamp_weight/samp_prob/mean(resamp_weight/samp_prob)

  MSE_SSL <- SemiSupervisedAccuracy(basis_labeled, basis_unlabeled,
                                    X_labeled, X_unlabaled, y, beta_SSL,
                                    beta_imp, weight = weight, threshold = NULL)

  OMR_SSL <- SemiSupervisedAccuracy(basis_labeled, basis_unlabeled,
                                    X_labeled, X_unlabaled, y, beta_SSL,
                                    beta_imp, weight = weight,
                                    threshold = threshold)

  return(list(mse_ssl = MSE_SSL, omr_ssl = OMR_SSL))
}
