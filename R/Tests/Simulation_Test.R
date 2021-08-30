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

# Cross-validated residuals.
num_folds <- 3
cv_residuals <- CrossValResids(basis_labeled, basis_unlabeled, X_labeled,
                               X_unlabeled, y, samp_prob, num_folds)

# Standard error estimates for supervised and SS estimates of regression parameter
beta_SL_se <- StdErrorEstimation(X_labeled, X_unlabeled, y,
                                 beta_sl, cv_residuals$beta_sl_cv)
beta_SSL_se <- StdErrorEstimation(X_labeled, X_unlabeled, y,
                                  beta_ssl, cv_residuals$beta_ssl_cv)

# Minimum Variance Estimator (here the component-wise optimal estimator)
# constant for stability due to high correlation between beta_SL and beta_SSL
eps.s <- (n.t*(se.beta.sl^2 + se.beta.ssl^2))/(2*n.t^0.6)
beta.ssl.w.ob <- min.var.est(beta.ssl, beta.sl, resids.gamma,
                            resids.beta.sl, Xt, Xv, eps.s);

# Our ensemble estimator (of the SL and SSL) for beta:
beta.ssl.w = beta.ssl.w.ob$beta

w.beta = beta.ssl.w.ob$weight
se.beta.ssl.w = beta.ssl.w.ob$se.est

# Estimated asymptotic standard error of the estimators (used for confidence interval construction):
se.beta.sl
se.beta.ssl
se.beta.ssl.w
se.beta.dr

################################################################################




# Run cross-validation.




