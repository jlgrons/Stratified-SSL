# Updated: 2021-04-22

#' Apparent estimates for brier score and overall misclassification rate.
#'
#' @param basis_labeled Basis matrix for labeled data set.
#' @param bais_unlabeled Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of weights.
#' @param weight Numeric vector of resampling weights.
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
  weight <- resamp_weight/samp.prob/mean(resamp_weight/samp.prob)
  X_all <- rbind(X_labeled, X_unlabeled)

  pred_prob_SSL_all <- Logit(cbind(1,X_all) %*% beta_ssl)
  classification_SSL_all <- I(pred_prob_SSL_all > c)

  # Note: Used to be a typo here for beta SL.
  pred_prob_SSL_labeled <- Logit(cbind(1,X_labeled) %*% beta_ssl)
  classification_SSL_labeled <- I(pred_prob_SSL_labeled > c)

  imps_labeled <- cbind(1, basis_labeled)%*%beta_imp

  refit_MSE <- glm(Yt~cbind(lp.t.ssl), offset = imps.t, family = 'binomial', weights = weight)$coeff
  imps_MSE <- g.logit(cbind(1,pred_prob_SSL_all)%*%refit.p1 +  cbind(1, basis.x)%*%gamma)
  MSE_SSL <- mean(imps.pe1 + (pred_prob_SSL_all - 2*imps.pe1)*pred_prob_SSL_all)

  refit_OMR <- glm(Yt~cbind(lp.t.ssl.ind), offset = imps.t, family = 'binomial', weights = weight)$coeff
  imps_OMR <- g.logit(cbind(1,classification_SSL_all)%*%refit.p2 +  cbind(1, basis.x)%*%gamma)
  OMR_SSL <- mean(imps.pe2 + (classification_SSL_all - 2*imps.pe2)*classification_SSL_all)

  return(list(mse_ssl = SSL_MSE, omr_sl = OMR_SSL))
}
