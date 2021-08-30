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
  resamp_weight <- sapply(1:num_resamples, function(kk) 4*rbeta(n_labeled,
                                                              0.5, 1.5))

  resids_beta_SL_weighted <- ((resamp_weight - 1) * resids_beta_SL)
  resids_beta_imp_weighted <- ((resamp_weight - 1) * resids_beta_imp)

  # Note: double check correctness of this.
  X_labeled_intercept <- cbind(1, X_labeled)
  IF_beta_imp <- inverse_information %*% t(
    X_labeled_intercept) %*% resids_beta_imp_weighted / n_labeled
  IF_beta_SL <- inverse_information %*% t(
    X_labeled_intercept) %*% resids_beta_SL_weighted / n_labeled


  beta_SSL_pert <- beta_SSL + diag(min_var_weight) %*% IF_beta_SL + diag(
    1-min_var_weight) %*% IF_beta_imp
  beta_SL_pert <- beta_SL + IF_beta_imp

  beta_imp_pert <- sapply(1:num_resamples, function(kk){
    RidgeRegression(basis_labeled, y,
                    weights = resamp_weight[,kk] / samp_prob/mean(
                      resamp_weight[,kk] / samp_prob),
             lambda = log(ncol(basis_labeled)) / (n_labeled^1.5))
  })

  perturbations_sl <- lapply(1:num_resamples, function(kk){
    SupervisedApparentAccuracy(X_labeled, y, beta_SL_pert[,kk], samp_prob,
                               resamp_weight = resamp_weight[ ,kk],
                               threshold = my_threshold)})

  perturbations_ssl <- lapply(1:num_resamples, function(kk){
    SemiSupervisedApparentAccuracy(basis_labeled, basis_unlabeled,
                                   X_labeled, X_unlabaled,
                                   y, beta_SSL_pert[,kk],  beta_imp_pert[, kk],
                                   samp_prob,
                                   resamp_weight = resamp_weight[ ,kk],
                                   threshold = my_threshold)})

  ssl_pert_mse <- sapply(perturbations_ssl, "[[", 1);
  ssl_pert_ae <- sapply(perturbations_ssl, "[[", 2);

  sl_pert_mse <- sapply(perturbations_sl, "[[", 1);
  sl_pert_ae <- sapply(perturbations_sl, "[[", 2);

  return(list(ssl_pert_mse = ssl_pert_mse,
              sl_pert_mse = sl_pert_mse,
              ssl_pert_ae = ssl_pert_ae,
              sl_pert_ae = sl_pert_ae,
              beta_imp_pert = beta_imp_pert,
              beta_SSL_pert = beta_SSL_pert,
              beta_SL_pert= beta_SL_pert))

}