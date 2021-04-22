# Updated: 2021-04-22

#' Apparent estimates for brier score (MSE) and misclassification rate (OMR).
#'
#' @param X_labeled Covariate matrix for labeled data set.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of weights.
#' @param resamp_weight Numeric vector of resampling weights.
#' @param threshold Threshold for overall misclassification rate.
#' @export
#' @return Supervised MSE and OMR.
#'

SemiSupervisedApparentAccuracy <- function(X_labeled, y, beta_SSL,
                                           beta_imp, samp_prob,
                                           resamp_weight = NULL,
                                           threshold = 0.5){
# computes the apparent estimates for model evaluation parameters
SupervisedApparentAccuracy <- function(basis_labeled, basis_unlabeled,
                                           X_labeled, X_unlabaled,
                                           y, beta_SL, beta_SSL, beta_imp,
                                           beta_DR
                                           Yt, beta.sl, beta.ssl, gamma, beta.dr,
                                           Xt, Xv, basis.x, samp.prob, V = NULL, c = 0.5){

  ## beta.sl: supervised beta
  ## beta.ssl: semi-supervised beta
  ## gamma: estimated parameter for imputation model
  ## Xt: covariates for labeled
  ## Xv: covariates for unlabeled
  ## Yt: labels
  ## basis.x: basis exp matrix
  ## samp.prob: sampling probability
  ## V: optional weights for resampling
  ## c: cutoff for the overall misclassification rate

  if(is.null(V)){V = rep(1, length(Yt))};
  weight = V/samp.prob/mean(V/samp.prob);
  dat.all = rbind(Xt, Xv);

  lp.v.ssl = g.logit(cbind(1,dat.all)%*%beta.ssl);
  lp.v.ssl.ind = I(lp.v.ssl > c);

  lp.t.ssl = g.logit(cbind(1,Xt)%*%beta.sl);
  lp.t.ssl.ind = I(lp.t.ssl > c);

  lp.t.dr = g.logit(cbind(1,Xt)%*%beta.dr);

  lp.t.dr.ind = I(lp.t.dr > c);

  mse.sl = mean((Yt - lp.t.ssl)^2 * weight)
  ae.sl = mean(abs(Yt - lp.t.ssl.ind) * weight)

  mse.naive = mean((Yt - lp.t.ssl)^2)
  ae.naive = mean(abs(Yt - lp.t.ssl.ind))

  mse.dr = mean((Yt - lp.t.dr)^2 * weight)
  ae.dr = mean(abs(Yt - lp.t.dr.ind) * weight)

  return(list(mse.ssl = mse.ssl, ae.ssl = ae.ssl,
              mse.sl = mse.sl, ae.sl = ae.sl, mse.dr = mse.dr,
              ae.dr = ae.dr, mse.naive = mse.naive, ae.naive = ae.naive))
}
