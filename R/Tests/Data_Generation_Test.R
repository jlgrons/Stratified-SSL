## Simulations for the Main Manuscript ##

## CC: Correct Prediction and Imputation Model

# OLD CODE
setwd('~/Desktop/Stratified-SSL/R/Original Code')
source('function.R')
set.seed(92047)
# Sample size of label data
n.t <- 400
# Matrix correlation
rho <- 0.4
# Sample size of unlabeled data
N <- 20000
# Number of covariates
p <- 10
# Number of stratum
strata.num <- 2
# Generate Data
old_data <- data_gen(n.t, N, p, rho, model.spec = 'CC', strata_num = strata.num)

# NEW CODE
setwd('~/Desktop/Stratified-SSL/R')
source('DataGeneration.R')
set.seed(92047)
n_lab <- 400
n_unlab <- 20000
p <- 10
rho <- 0.4
num_strata <- 2
new_data <- DataGeneration(n_lab, n_unlab, p, rho, signal = c(1, 1, 0.5, 0.5),
                           model_specification = 'outcome_correct_imputation_correct',
                           num_strata = num_strata)

# Test to see if they output same data
all.equal(old_data$b0, new_data$signal)

all.equal(old_data$X0, new_data$covariates)
all.equal(old_data$Xt, new_data$covariates_lab)
all.equal(old_data$Xv, new_data$covariates_unlab)
all.equal(old_data$Xr, new_data$covariates_random_samp)
all.equal(old_data$Xvr, new_data$covariates_unlab_random_samp)

all.equal(old_data$Yt, new_data$Y_lab)
all.equal(old_data$Y, new_data$Y)
all.equal(old_data$Yr, new_data$Y_random_samp)

all.equal(old_data$St, new_data$S_lab)
all.equal(old_data$Sv, new_data$S_unlab)
all.equal(old_data$Sr, new_data$S_random_samp)
all.equal(old_data$S, new_data$strat_var)
all.equal(old_data$samp.prob, new_data$samp_prob)

## IC


## II


## Simulations for the Supplementary Materials ##

## GM


## IC1


## II1

