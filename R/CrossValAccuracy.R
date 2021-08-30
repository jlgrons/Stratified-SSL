# Updated: 2021-08-30

#' CV estimates for brier score (MSE) and misclassification rate (OMR).
#'
#' @param basis_labeled Basis matrix for labeled data set.
#' @param bais_unlabeled Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of weights.
#' @param min_var_weight Numeric vector of minimum variance weights.
#' @param num_folds Scalar indicating number of folds for CV.
#' @param num_folds Scalar indicating number of repitions for CV.
#' @param threshold Threshold for overall misclassification rate.
#' @param lambda0 Initial lambda for imputation model.
#' @export
#' @return CV semi-supervised and supervised MSE and OMR.
#'

# computes the CV estimates for model evaluation parameters
CrossValAccuaracy <- function(basis_labeled, basis_unlabeled,
                              X_labeled, X_unlabeled, y, samp_prob,
                              min_var_weight, num_folds = 3, reps = 10,
                              theshold = 0.5, lambda0 = NULL){

  pp <- ncol(basis_labeled)
  n_labeled <- nrow(X_labeled)
  data.all <- rbind(X_labeled,X_unlabeled)

  mse_cv_ssl <- matrix(NA, nrow = num_folds, ncol = reps)
  ae_cv_ssl <- matrix(NA, nrow = num_folds, ncol = reps)

  mse_cv_sl <- matrix(NA, nrow = num_folds, ncol = reps)
  ae_cv_sl <- matrix(NA, nrow = num_folds, ncol = reps)

  mse_cv_naive <- matrix(NA, nrow = num_folds, ncol = reps)
  ae_cv_naive <- matrix(NA, nrow = num_folds, ncol = reps)

  mse_cv_dr <- matrix(NA, nrow = num_folds, ncol = reps)
  ae_cv_dr <- matrix(NA, nrow = num_folds, ncol = reps)


  if(is.null(lambda0)){lambda0 = log(pp)/(floor((num_folds-1)*n_labeled/num_folds))^1.5}

  for(j in 1:reps){
    set.seed(j)
    ind.cv = split(1:n_labeled, sample(reps(1:num_folds, floor(n_labeled/num_folds))))

    for(i in 1:num_folds){

      inds.v = as.vector(unlist(ind.cv[i]))
      inds.t = as.vector(setdiff(as.vector(unlist(ind.cv)), inds.v))

      wg.v = samp.prob[inds.v] wg.t = samp.prob[inds.t]

      Yt.v = Yt[inds.v] X_labeled.v = X_labeled[inds.v,] basis.v = basis.x[inds.v, ]
      Yt.t = Yt[inds.t] X_labeled.t = X_labeled[inds.t,] basis.t = basis.x[inds.t, ]

      beta.tmp = glm.fit.SS(rbind(basis.t, basis.x[-ind.lab, ]),
                            X_labeled.t, X_unlabeled, Yt.t, wg.t, lambda0)

      beta.ssl.t =  beta.tmp$beta.ssl
      beta.sl.t = beta.tmp$beta.sl
      gamma.t =  beta.tmp$gamma
      beta.ssl.w.t = W*beta.ssl.t + (1-W)*beta.sl.t
      beta.dr.t = beta.tmp$beta.dr

      lp.v = g.logit(cbind(1,X_labeled.v)%*%beta.ssl.w.t)
      lp.v.ind = ifelse(I(lp.v > c), 1, 0)

      imps.v = cbind(1, basis.v)%*%gamma.t
      refit.pv.1 <- glm(Yt.v~cbind(lp.v), offset = imps.v,
                        weights = 1/wg.v/mean(1/wg.v), family = 'binomial')$coeff
      refit.pv.2 <-glm(Yt.v~cbind(lp.v.ind), offset = imps.v,
                       weights = 1/wg.v/mean(1/wg.v), family = 'binomial')$coeff

      lp.u = g.logit(cbind(1,data.all)%*%beta.ssl.w.t)
      lp.u.ind = I(lp.u > c)

      imps.pe.1 = g.logit(cbind(1,lp.u)%*% refit.pv.1 + cbind(1, basis.x)%*%gamma.t)
      imps.pe.2 = g.logit(cbind(1,lp.u.ind)%*% refit.pv.2 + cbind(1, basis.x)%*%gamma.t)

      mse_cv_ssl[i,j] = mean(imps.pe.1 + (lp.u - 2*imps.pe.1)*lp.u)
      ae_cv_ssl[i,j] = mean(imps.pe.2 + (lp.u.ind - 2*imps.pe.2)*lp.u.ind)

      lp.v.sl = g.logit(cbind(1,X_labeled.v) %*% beta.sl.t)
      lp.v.sl.ind = I(lp.v.sl > c)

      mse_cv_sl[i,j] = mean((Yt.v  - lp.v.sl)^2 * 1/wg.v)/mean(1/wg.v)
      ae_cv_sl[i,j] = mean(abs(Yt.v - lp.v.sl.ind) * 1/wg.v)/mean(1/wg.v)

      mse_cv_naive[i,j] = mean((Yt.v  - lp.v.sl)^2)
      ae_cv_naive[i,j] = mean(abs(Yt.v - lp.v.sl.ind))

      lp.v.dr = g.logit(cbind(1, X_labeled.v) %*% beta.dr.t)
      lp.v.dr.ind = I(lp.v.dr > c)

      mse_cv_dr[i,j] = mean((Yt.v  - lp.v.dr)^2 * 1 / wg.v) / mean(1 / wg.v)
      ae_cv_dr[i,j] = mean(abs(Yt.v - lp.v.dr.ind) * 1 / wg.v) / mean(1 / wg.v)

    }
  }

  return(list(mse.ssl = mean(mse_cv_ssl, na.rm = T),
              ae.ssl = mean(ae_cv_ssl, na.rm = T),
              mse.sl = mean(mse_cv_sl, na.rm = T),
              ae.sl = mean(ae_cv_sl, na.rm = T),
              mse.dr = mean(mse_cv_dr, na.rm = T),
              ae.dr = mean(ae_cv_dr, na.rm = T),
              mse.naive = mean(mse_cv_naive, na.rm = T),
              ae.naive = mean(ae_cv_naive, na.rm = T)))

}
