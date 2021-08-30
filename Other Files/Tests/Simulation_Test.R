################################################################################
# Simulation setting test
library('stratifiedSSL')
set.seed(92047)

n_lab <- 400
n_unlab <- 20000
p <- 10
rho <- 0.4
num_strata <- 2

################################################################################
# Generate data.
new_data <- DataGeneration(n_lab, n_unlab, p, rho, signal = c(1, 1, 0.5, 0.5),
                           model_specification =
                             'outcome_incorrect_imputation_correct',
                           num_strata = num_strata)

# Format data.
X_labeled <- new_data$covariates_lab
X_unlabeled <- new_data$covariates_unlab
S_unlabeled <- new_data$S_unlab
S_labeled <- new_data$S_lab
y <- new_data$Y_lab
samp_prob <- new_data$samp_prob

################################################################################
# Get basis expansion.
my_basis <- NaturalSplineBasis(rbind(X_labeled, X_unlabeled),
                               c(S_labeled, S_unlabeled),
                               num_knots = 3)

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

# Standard error estimates for supervised and SS estimates.
beta_sl_se_obj <- StdErrorEstimation(X_labeled, X_unlabeled, y,
                                 beta_sl, resids_sl)
beta_ssl_se_obj <- StdErrorEstimation(X_labeled, X_unlabeled, y,
                                  beta_ssl, resids_ssl)
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

################################################################################
# Apparent accuracy estimates.
my_threshold <- 0.5
acc_sl <- SupervisedApparentAccuracy(X_labeled, y, beta_sl, samp_prob,
                                       resamp_weight = NULL,
                                       threshold = my_threshold)

acc_ssl <- SemiSupervisedApparentAccuracy(basis_labeled, basis_unlabeled,
                                          X_labeled, X_unlabaled,
                                          y, beta_ssl, beta_imp, samp_prob,
                                          resamp_weight = NULL,
                                          threshold = my_threshold)

################################################################################
# Cross-validated accuracy estimates.
reps <- 1
# Note: Function still needs to be formatted.
acc_cv <- CrossValAccuracy(basis_labeled, basis_unlabeled,
                           X_labeled, X_unlabeled, y, samp_prob,
                           beta_minvar$min_var_weight[,1],
                           num_folds = num_folds, reps = reps,
                           theshold = my_threshold, lambda0 = NULL)

# Ensemble of apparent and CV estimator.
#w <- K.fold / (2 * K.fold - 1)
#w * ap.ae + (1 - w) * cv.ae
#w * ap.mse + (1 - w) * cv.mse

################################################################################
# Perturbation for standard error estimates for accuracy estimates.
num_perts <- 2
acc_pert <- AccuracyStdErrorEstimation(basis_labeled, basis_unlabeled,
                                       X_labeled, X_unlabeled, y,
                                       samp_prob, beta_minvar$min_var_weight[,1],
                                       beta_sl, beta_ssl,
                                       cv_residuals$resids_beta_sl,
                                       cv_residuals$resids_beta_imp,
                                       beta_ssl_se$inverse_information,
                                       num_resamples = num_perts,
                                       threshold = my_threshold)

