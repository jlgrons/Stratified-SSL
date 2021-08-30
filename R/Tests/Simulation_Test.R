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

################################################################################




# Run cross-validation.




