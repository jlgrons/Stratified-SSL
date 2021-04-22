# Updated: 2021-04-19

#' Component-wise minimum variance semi-supervised regression.
#'
#' @param basis_labeled Basis matrix for labeled data set.
#' @param bais_unlabeled Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of weights.
#' @param min_var_weight Minimum variance weight for semi-supervised estimate.
#' @param beta_SL Supervised regression coefficient vector.
#' @param beta_SSL Semi-supervised regression coefficient vector.
#' @param resids_beta_SL Residuals from the supervised regression model.
#' @param resids_beta_SSL Residuals from the semi-supervised regression model.
#' @param resids_beta_imp Residuals from the imputation model.
#' @param inverse_information Inverse information  matric.
#' @param num_resamples Number of resamples.
#' @param threshold Threshold for over misclassification rate.
#' @export
#' @return Pertrubed estimates.
#'

AccuracyStdErrorEstimation <- function(basis_labeled, basis_unlabeled,
                                       X_labeled, X_unlabeled, y,
                                       samp_prob, min_var_weight,
                                       beta_SL, beta_SSL, resids_beta_SL,
                                       resids_beta_imp, inverse_information,
                                       num_resamples = 500, threshold = 0.5){

  n_labeled <- length(y)
  resamp_weight <- sapply(num_resamples, function(kk) 4*rbeta(n_labeled, 0.5, 1.5))

  resids_beta_SL_weighted <- ((resamp_weight - 1) * resids_beta_SL)
  resids_beta_imp_weighted <- ((resamp_weight - 1) * resids_beta_imp)

  # Note: double check correctness of this.
  X_labeled_intercept <- cbind(1, X_labeled)
  IF_beta_imp <- inverse_information %*% t(
    X_labeled_intercept) %*% resids_beta_imp_weighted / n_labeled
  IF_beta_SL <- inverse_information %*% t(
    X_labeled_intercept) %*% resids_beta_SL_weighted / n_labeled


  beta_SSL_pert <- beta_SSL + diag(min_var_weight[,1]) %*% IF_beta_SL + diag(
    min_var_weight[,2]) %*% IF_beta_imp
  beta_SL_pert <- beta_SL + IF_beta_imp

  beta_imp_pert <- sapply(num_resamples, function(kk){
    my.ridge(basis.x[ind.lab, ], Yt, weights = resamp_weight[,kk]/samp.prob/mean(resamp_weight[,kk]/samp.prob),
             lambda = log(ncol(basis_labeled))/(n_labeled^1.5))
  })

  perts <- lapply(num_resamples, function(kk){model.eval.ap(Yt, beta.sl.p[,kk], beta_SSL_pert[,kk], beta_imp_pert[,kk],
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
              beta_imp_pert = beta_imp_pert, beta_SSL_pert = beta_SSL_pert,
              beta_SL_pert= beta.sl.p, beta.dr.p = beta.dr.p))

}
