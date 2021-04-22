# Updated: 2021-04-22

#' Apparent estimates for brier score and overall misclassification rate.
#'
#' @param basis_labeled Basis matrix for labeled data set.
#' @param bais_unlabeled Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of weights.
#' @param lambda Penalization parameter for initial ridge estimator.
#' @export
#' @return Vector containing regression coefficients.
#'

# computes the apparent estimates for model evaluation parameters
SemiSupervisedApparentAccuracy <- function(basis_labeled, basis_unlabeled,
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

  n.t = length(Yt); ind.l = (1:n.t);
  imps.t = cbind(1, basis.x[ind.l,])%*%gamma

  refit.p1 <- glm(Yt~cbind(lp.t.ssl), offset = imps.t, family = 'binomial', weights = weight)$coeff;
  refit.p2 <- glm(Yt~cbind(lp.t.ssl.ind), offset = imps.t, family = 'binomial', weights = weight)$coeff;
  imps.pe1 = g.logit(cbind(1,lp.v.ssl)%*%refit.p1 +  cbind(1, basis.x)%*%gamma);
  imps.pe2 = g.logit(cbind(1,lp.v.ssl.ind)%*%refit.p2 +  cbind(1, basis.x)%*%gamma);

  mse.ssl = mean(imps.pe1 + (lp.v.ssl - 2*imps.pe1)*lp.v.ssl)
  ae.ssl = mean(imps.pe2 + (lp.v.ssl.ind - 2*imps.pe2)*lp.v.ssl.ind)

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
