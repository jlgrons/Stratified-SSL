# Updated: 2021-04-22

#' Semi-supervised accuracy estimation.
#'
#' @param basis_labeled Basis matrix for labeled data set.
#' @param bais_unlabeled Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param beta_SSL Numeric vector of regression coefficients.
#' @param beta_imp Numeric vector of regression coefficients for imputation.
#' @param ind_est Numeric vector to indicate inds for refitting model.
#' @param weight Numeric vector of observation weights.
#' @param threshold Threshold for overall misclassification rate.
#' @export
#' @return Semi-supervised accuracy.
#'

SemiSupervisedAccuracy <- function(basis_labeled, basis_unlabeled,
                                   X_labeled, X_unlabaled, y, beta_SSL,
                                   beta_imp, ind_est,
                                   weight = NULL, threshold = NULL){

  if(is.null(weight)){weight <- rep(1, length(y))}

  X_all <- rbind(X_labeled, X_unlabeled)
  basis_all <- rbind(basis_labeled, basis_unlabeled)

  if(is.null(threshold)){
    pred_prob_SSL_all <- Expit(cbind(1, X_all) %*% beta_SSL)
    pred_prob_SSL_labeled <- Expit(cbind(1, X_labeled) %*% beta_SSL)
  }else{
    pred_prob_SSL_all <- I(Expit(cbind(1, X_all) %*% beta_SSL) > threshold)
    pred_prob_SSL_labeled <- I(Expit(cbind(1, X_labeled) %*% beta_SSL) > threshold)
  }

  imps_labeled <- cbind(1, basis_labeled) %*% beta_imp
  imps_all <- cbind(1, basis_all) %*% beta_imp

  # Fit only using specified indices for CV.
  beta_refit <- glm(y[ind_est] ~ cbind(pred_prob_SSL_labeled)[ind_est],
                    offset = imps_labeled[ind_est],
                    family = 'binomial', weights = weight)$coeff
  
  # Project to the full data.
  imps_refit <- Expit(cbind(1, pred_prob_SSL_all) %*% beta_refit +
                          imps_all)

  accuracy_SSL <- mean(imps_refit +
                         (pred_prob_SSL_all - 2*imps_refit)*pred_prob_SSL_all)

  return(accuracy_SSL)
}
