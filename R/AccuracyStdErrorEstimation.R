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

  n_labeled <- length(y)
  resamp_weight <- sapply(num_resamples, function(kk) 4*rbeta(n_labeled, 0.5, 1.5))

  resids_beta_sl_weighted <- ((resamp_weight - 1) * resids_beta_sl)
  resids_beta_imp_weighted <- ((resamp_weight - 1) * resids_beta_imp)

  # Note: double check correctness of this.
  X_labeled_intercept <- cbind(1, X_labeled)
  IF_beta_imp <- inverse_information %*% t(
    X_labeled_intercept) %*% resids_beta_imp_weighted / n_labeled
  IF_beta_sl <- inverse_information %*% t(
    X_labeled_intercept) %*% resids_beta_sl_weighted / n_labeled


  beta_ssl_pert <- beta_ssl + diag(min_var_weight[,1]) %*% IF_beta_sl + diag(
    min_var_weight[,2]) %*% IF_beta_imp
  beta_sl_pert <- beta_sl + IF_beta_imp

  beta_imp_pert <- sapply(num_resamples, function(kk){
    my.ridge(basis.x[ind.lab, ], Yt, weights = resamp_weight[,kk]/samp.prob/mean(resamp_weight[,kk]/samp.prob),
             lambda = log(ncol(basis_labeled))/(n_labeled^1.5))
  })

  perts <- lapply(num_resamples, function(kk){model.eval.ap(Yt, beta.sl.p[,kk], beta_ssl_pert[,kk], beta_imp_pert[,kk],
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
              beta_imp_pert = beta_imp_pert, beta_ssl_pert = beta_ssl_pert,
              beta_sl_pert= beta.sl.p, beta.dr.p = beta.dr.p))

}
