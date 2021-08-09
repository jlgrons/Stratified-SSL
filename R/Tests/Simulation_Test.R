# Simulation setting test


set.seed(92047)

n_lab <- 400
n_unlab <- 20000
p <- 10
rho <- 0.4
num_strata <- 2


new_data <- DataGeneration(n_lab, n_unlab, p, rho, signal = c(1, 1, 0.5, 0.5),
                           model_specification = 'outcome_incorrect_imputation_correct',
                           num_strata = num_strata)

X_labeled <- new_data$covariates_lab
X_unlabeled <- new_data$covariates_unlab



regression <- SemiSupervisedRegression(basis_labeled,
                                       basis_unlabeled,
                                       X_labeled,
                                       X_unlabeled,
                                       y,
                                       samp_prob,
                                       lambda = NULL)
