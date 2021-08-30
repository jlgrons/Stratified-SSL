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
                           model_specification = 'outcome_incorrect_imputation_correct',
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
# Beta from imputation.
beta_imp <- regression_result$beta_imp

# Cross-validated residuals.
num_folds <- 3
cv_residuals <- CrossValResids(basis_labeled, basis_unlabeled, X_labeled,
                               X_unlabeled, y, samp_prob, num_folds)

# Standard error estimates for supervised and SS estimates.
beta_sl_se <- StdErrorEstimation(X_labeled, X_unlabeled, y,
                                 beta_sl, cv_residuals$beta_sl_cv)
beta_ssl_se <- StdErrorEstimation(X_labeled, X_unlabeled, y,
                                  beta_ssl, cv_residuals$beta_ssl_cv)

# Minimum Variance Estimator (here the component-wise optimal estimator).
my_epsilon <- (n_labeled*(beta_sl_se$std_error^2 +
                            beta_ssl_se$std_error^2))/(2*n_labeled^0.6)
beta_minvar <- SemiSupervisedMinVarRegression(beta_ssl, beta_sl,
                                              beta_ssl_se$std_error,
                                              beta_sl_se$std_error,
                                              cv_residuals$resids_beta_imp,
                                              cv_residuals$resids_beta_sl,
                                              X_labeled,
                                              X_unlabeled, epsilon = my_epsilon)

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

# Ensemble of apparent and CV estimator (Section 4)
w <- K.fold / (2 * K.fold - 1)
w * ap.ae + (1 - w) * cv.ae
w * ap.mse + (1 - w) * cv.mse
################################################################################
# Perturbation for standard error estimates for accuracy estimates.


