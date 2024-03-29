# This file will run the simulations presented in the manuscript.

################################################################################
# Load the R package.
library('stratifiedSSL')

# Set seed for replication.
set.seed(92047)

# Set desired parameters.
n_lab <- 400
n_unlab <- 20000
p <- 10
rho <- 0.4
num_strata <- 2

################################################################################
# Generate data.

# Three 'model_specification' settings from the main text:
# outcome_correct_imputation_correct
# outcome_incorrect_imputation_correct
# outcome_incorrect_imputation_incorrect

# Three 'model_specification' settings from the supplement:
# outcome_incorrect_imputation_correct_supp
# outcome_incorrect_imputation_incorrect_supp
# gaussian_mixture

# Note: The example file for the gaussian mixture setting will be added.  In the 
# interim, you can run the supp_GM.R file in Other Files > Original Code.

new_data <- DataGeneration(n_lab, n_unlab, p, rho, signal = c(1, 1, 0.5, 0.5),
                           model_specification =
                             'outcome_correct_imputation_correct',
                           num_strata = num_strata)

# Format data.
X_labeled <- new_data$covariates_lab
X_unlabeled <- new_data$covariates_unlab
S_unlabeled <- new_data$S_unlab
S_labeled <- new_data$S_lab
y <- new_data$Y_lab
samp_prob <- new_data$samp_prob

################################################################################
# Get basis expansion. Uses the basis specification from the paper.

# Three 'model_specification' settings from the main text:
# outcome_correct_imputation_correct: basis_type <- 'NS_basis'
# outcome_incorrect_imputation_correct: basis_type <- 'interact'
# outcome_incorrect_imputation_incorrect: basis_type <- 'NS_basis'

# Three 'model_specification' settings from the supplement:
# outcome_incorrect_imputation_correct_supp: basis_type <- 'IC1'
# outcome_incorrect_imputation_incorrect_supp: basis_type <- 'II1'

num_knots <- 3
basis_type <- 'NS_basis'
if(basis_type == 'NS_basis'){
  
  my_basis <- NaturalSplineBasis(rbind(X_labeled, X_unlabeled),
                                 c(S_labeled, S_unlabeled),
                                 num_knots = num_knots)
  
}else{
  
  my_basis <- AlternativeBasis(rbind(X_labeled, X_unlabeled),
                               c(S_labeled, S_unlabeled),
                               num_knots = num_knots, basis_type = basis_type) 
  
}

basis_labeled <- my_basis[1:n_lab, ]
basis_unlabeled <- my_basis[(n_lab+1):nrow(my_basis), ]

################################################################################
# Fit the regression model.
regression_result <- SemiSupervisedRegression(basis_labeled,
                                              basis_unlabeled,
                                              X_labeled,
                                              X_unlabeled,
                                              y,
                                              samp_prob,
                                              lambda = NULL)

# Supervised beta.
beta_sl <- regression_result$beta_SL
# Semi-supervised beta.
beta_ssl <- regression_result$beta_SSL
# Naive beta with no sampling probability.
beta_naive <- regression_result$beta_SL_unweighted
# Gamma from imputation.
gamma <- regression_result$gamma
# Beta from density ratio method.
beta_dr <- regression_result$beta_DR
# DR projection.
proj_dr <- regression_result$proj_DR

################################################################################
# Cross-validated residuals.
num_folds <- 3
cv_residuals <- CrossValResids(basis_labeled, basis_unlabeled, X_labeled,
                               X_unlabeled, y, samp_prob, num_folds)

# Save the CV residuals.
resids_sl <- cv_residuals$resids_beta_sl
resids_ssl <- cv_residuals$resids_beta_ssl
resids_dr <- cv_residuals$resids_beta_dr
resids_gamma <- cv_residuals$resids_gamma

# Save all the resids.
resids_all <- cbind(resids_ssl, resids_sl, resids_dr, resids_gamma)

# Standard error estimates for supervised and SS estimates.
beta_sl_se_obj <- StdErrorEstimation(X_labeled, X_unlabeled, y,
                                     beta_sl, resids_sl)
beta_ssl_se_obj <- StdErrorEstimation(X_labeled, X_unlabeled, y,
                                      beta_ssl, resids_gamma)
beta_dr_se_obj <- StdErrorEstimationDR(X_labeled, X_unlabeled, y,
                                       beta_dr, resids_dr, proj_dr)

# Supervised SE and inverse info.
beta_sl_se <- beta_sl_se_obj$std_error
beta_sl_inv_info <- beta_sl_se_obj$inverse_information
# Semi-supervised SE and inverse info.
beta_ssl_se <- beta_ssl_se_obj$std_error
beta_ssl_inv_info <- beta_ssl_se_obj$inverse_information
# Semi-supervised SE and inverse info.
beta_dr_se <- beta_dr_se_obj

# Minimum Variance Estimator (here the component-wise optimal estimator).
beta_minvar <- SemiSupervisedMinVarRegression(beta_ssl, beta_sl,
                                              beta_ssl_se,
                                              beta_sl_se,
                                              resids_gamma,
                                              resids_sl,
                                              X_labeled,
                                              X_unlabeled,
                                              epsilon = NULL)

# Save result.
beta_mv <- beta_minvar$beta_SSL_min_var
mv_weight <- beta_minvar$min_var_weight
beta_mv_se <- beta_minvar$se_beta_weight

# Save all the standard errors.
se_all <- cbind(beta_ssl_se,  beta_mv_se, beta_sl_se, beta_dr_se)

# Save all the betas.
beta_all <- cbind(beta_ssl, beta_mv, beta_sl, beta_dr, beta_naive)
################################################################################
# Apparent accuracy estimates.
my_threshold <- 0.5

acc_sl <- SupervisedApparentAccuracy(X_labeled, y, beta_sl, samp_prob,
                                     resamp_weight = NULL,
                                     threshold = my_threshold)

acc_ssl <- SemiSupervisedApparentAccuracy(basis_labeled, basis_unlabeled,
                                          X_labeled, X_unlabaled,
                                          y, beta_mv, gamma, samp_prob,
                                          resamp_weight = NULL,
                                          threshold = my_threshold)

acc_dr <- SupervisedApparentAccuracy(X_labeled, y, beta_dr, samp_prob,
                                     resamp_weight = NULL,
                                     threshold = my_threshold)

acc_naive <- SupervisedApparentAccuracy(X_labeled, y, beta_sl,
                                        rep(1, length(y)),
                                        resamp_weight = NULL,
                                        threshold = my_threshold)

acc_ap_omr <- c(acc_ssl$omr_ssl, acc_sl$omr_sl, acc_dr$omr_sl, acc_naive$omr_sl)
acc_ap_mse <- c(acc_ssl$mse_ssl, acc_sl$mse_sl, acc_dr$mse_sl, acc_naive$mse_sl)

acc_ap_omr
acc_ap_mse

################################################################################
# Cross-validated accuracy estimates.
reps <- 2 
# reps is set to be low to run quickly as an example. We recommend 10-20.

acc_cv <- CrossValAccuracy(basis_labeled, basis_unlabeled,
                           X_labeled, X_unlabeled, y, samp_prob,
                           mv_weight[,1],
                           num_folds, reps,
                           my_threshold, lambda0 = NULL)

acc_cv_omr <- c(acc_cv$omr_ssl, acc_cv$omr_sl, acc_cv$omr_dr, acc_cv$omr_naive)
acc_cv_mse <- c(acc_cv$mse_ssl, acc_cv$mse_sl, acc_cv$mse_dr, acc_cv$mse_naive)

acc_cv_omr
acc_cv_mse

################################################################################
# Ensemble of apparent and CV estimator.
cv_weight <- num_folds / (2 * num_folds - 1)
acc_en_omr <- cv_weight * acc_ap_omr + ((1-cv_weight) * acc_cv_omr)
acc_en_mse <- cv_weight * acc_ap_mse + ((1-cv_weight) * acc_cv_mse)

acc_en_omr
acc_en_mse

################################################################################
# Perturbation for standard error estimates for accuracy estimates.
num_perts <- 2
# num_perts is set to be low to run quickly as an example. We recommend 500.
acc_pert <- AccuracyStdErrorEstimation(basis_labeled, basis_unlabeled,
                                       X_labeled, X_unlabeled, y,
                                       samp_prob,
                                       mv_weight[,1],
                                       beta_sl,
                                       beta_mv,
                                       beta_dr,
                                       resids_sl,
                                       resids_gamma,
                                       resids_dr,
                                       proj_dr,
                                       beta_ssl_inv_info,
                                       num_resamples = num_perts,
                                       threshold = my_threshold)

# Save results.
ssl_pert_mse <- acc_pert$ssl_pert_mse
sl_pert_mse <- acc_pert$sl_pert_mse
dr_pert_mse <- acc_pert$dr_pert_mse

ssl_pert_omr <- acc_pert$ssl_pert_omr
sl_pert_omr <- acc_pert$sl_pert_omr
dr_pert_omr <- acc_pert$dr_pert_omr

# Save all the results.
pert_mse_all <- cbind(ssl_pert_mse, sl_pert_mse, dr_pert_mse)
pert_omr_all <- cbind(ssl_pert_omr, sl_pert_omr, dr_pert_omr)

