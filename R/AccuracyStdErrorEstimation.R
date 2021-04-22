# Updated: 2021-04-19

#' Component-wise minimum variance semi-supervised regression.
#'
#' @param beta_SSL Semi-supervised regression coefficient vector.
#' @param beta_SL Supervised regression coefficient vector.
#' @param resids_imp Residuals from the imputation model.
#' @param resids_beta_SL Residuals from the supervised regression model.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param epsilon Small offset to help with numerical stability.
#' @export
#' @return List with parameter and SE est as well as estimated min_var_weights.
#'
#

AccuracyStdErrorEstimation <- function(basis_labeled, basis_unlabeled,
                                       X_labeled, X_unlabeled, y,
                                       samp_prob, min_var_weight,
                                       beta_sl, beta_ssl, resids_beta_sl,
                                       resids_beta_imp, inverse_information,
                                       num_resamples = 500, threshold = 0.5){

  resamp_weight <- sapply(num_resamples, function(kk) 4*rbeta(n.t, 0.5, 1.5))

  resids.beta.sl.p <- ((resamp_weight-1)*resids.beta.sl)
  resids.gamma.p <- ((resamp_weight-1)*resids.gamma)
  resids.beta.dr.p <- (resamp_weight-1)*resids.beta.dr
  proj.dr.p <- t(proj.dr) %*% (resamp_weight-1)

  Xt.1 <- cbind(1, Xt)
  T_1.p <- A %*% t(Xt.1)%*%resids.gamma.p/n.t;
  T_2.p <- A %*% t(Xt.1)%*%resids.beta.sl.p/n.t;

  T_1.dr.p <- A %*% t(Xt.1) %*% resids.beta.dr.p / n.t
  T_2.dr.p <- A %*% proj.dr.p / n.t

  beta.ssl.w.p <- beta.ssl.w + diag(W[,1]) %*% T_1.p + diag(W[,2]) %*% T_2.p;
  beta.sl.p <- beta.sl + T_2.p
  beta.dr.p <- beta.dr + T_1.dr.p - T_2.dr.p

  ind.lab <- 1:n.t;

  gamma.p <- sapply(num_resamples, function(kk){
    my.ridge(basis.x[ind.lab, ], Yt, weights = resamp_weight[,kk]/samp.prob/mean(resamp_weight[,kk]/samp.prob),
             lambda = log(ncol(basis.x))/(n.t^1.5))
  })

  perts <- lapply(num_resamples, function(kk){model.eval.ap(Yt, beta.sl.p[,kk], beta.ssl.w.p[,kk], gamma.p[,kk],
                                                 beta.dr.p[,kk], Xt, Xv, basis.x, samp.prob, V = resamp_weight[,kk])})

  ssl.pert.mse <- sapply(perts, "[[", 1);
  ssl.pert.ae <- sapply(perts, "[[", 2);

  sl.pert.mse <- sapply(perts, "[[", 3);
  sl.pert.ae <- sapply(perts, "[[", 4);

  dr.pert.mse <- sapply(perts, "[[", 5);
  dr.pert.ae <- sapply(perts, "[[", 6);

  return(list(ssl.pert.mse = ssl.pert.mse, sl.pert.mse = sl.pert.mse,
              dr.pert.mse = dr.pert.mse, ssl.pert.ae = ssl.pert.ae,
              sl.pert.ae = sl.pert.ae, dr.pert.ae = dr.pert.ae,
              gamma.p = gamma.p, beta.ssl.w.p = beta.ssl.w.p,
              beta.sl.p = beta.sl.p, beta.dr.p = beta.dr.p))

}
