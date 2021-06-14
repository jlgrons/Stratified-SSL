source('Testing_Functions.R')

## Simulations for the Main Manuscript ##

#################################################
## CC: Correct Prediction and Imputation Model ##
#################################################

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

# RUN TEST
test_result <- run_data_gen_test(old_data, new_data)
# Check if all tests passed
sum(unlist(test_result)) == length(test_result)

###########################################################
## IC: Incorrect Prediction and Correct Imputation Model ##
###########################################################

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
old_data <- data_gen(n.t, N, p, rho, model.spec = 'IC', strata_num = strata.num)

# NEW CODE
setwd('~/Desktop/Stratified-SSL/R')
source('DataGeneration.R')
source('HelperFunctions.R')
set.seed(92047)
n_lab <- 400
n_unlab <- 20000
p <- 10
rho <- 0.4
num_strata <- 2
new_data <- DataGeneration(n_lab, n_unlab, p, rho, signal = c(1, 1, 0.5, 0.5),
                           model_specification = 'outcome_incorrect_imputation_correct',
                           num_strata = num_strata)

# RUN TEST
test_result <- run_data_gen_test(old_data, new_data)
# Check if all tests passed
sum(unlist(test_result)) == length(test_result)


###################################################
## II: Incorrect Prediction and Imputation Model ##
###################################################

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
old_data <- data_gen(n.t, N, p, rho, model.spec = 'II', strata_num = strata.num)

# NEW CODE
setwd('~/Desktop/Stratified-SSL/R')
source('DataGeneration.R')
source('HelperFunctions.R')
set.seed(92047)
n_lab <- 400
n_unlab <- 20000
p <- 10
rho <- 0.4
num_strata <- 2
new_data <- DataGeneration(n_lab, n_unlab, p, rho, signal = c(1, 1, 0.5, 0.5),
                           model_specification = 'outcome_incorrect_imputation_incorrect',
                           num_strata = num_strata)

# RUN TEST
test_result <- run_data_gen_test(old_data, new_data)
# Check if all tests passed
sum(unlist(test_result)) == length(test_result)


# __________________________________________________________________________ #

## Simulations for the Supplementary Materials ##

###########################################################
## IC: Incorrect Prediction and Correct Imputation Model ##
###########################################################

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
old_data <- data_gen(n.t, N, p, rho, model.spec = 'IC1', strata_num = strata.num)

# NEW CODE
setwd('~/Desktop/Stratified-SSL/R')
source('DataGeneration.R')
source('HelperFunctions.R')
set.seed(92047)
n_lab <- 400
n_unlab <- 20000
p <- 10
rho <- 0.4
num_strata <- 2
new_data <- DataGeneration(n_lab, n_unlab, p, rho, signal = c(1, 1, 0.5, 0.5),
                           model_specification = 'outcome_incorrect_imputation_correct_supp',
                           num_strata = num_strata)

# RUN TEST
test_result <- run_data_gen_test(old_data, new_data)
# Check if all tests passed
sum(unlist(test_result)) == length(test_result)

###################################################
## II: Incorrect Prediction and Imputation Model ##
###################################################

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
old_data <- data_gen(n.t, N, p, rho, model.spec = 'II1', strata_num = strata.num)

# NEW CODE
setwd('~/Desktop/Stratified-SSL/R')
source('DataGeneration.R')
source('HelperFunctions.R')
set.seed(92047)
n_lab <- 400
n_unlab <- 20000
p <- 10
rho <- 0.4
num_strata <- 2
new_data <- DataGeneration(n_lab, n_unlab, p, rho, signal = c(1, 1, 0.5, 0.5),
                           model_specification = 'outcome_incorrect_imputation_incorrect_supp',
                           num_strata = num_strata)

# RUN TEST
test_result <- run_data_gen_test(old_data, new_data)
# Check if all tests passed
sum(unlist(test_result)) == length(test_result)


#############################
## Gaussian Mixture Model ##
############################

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
old_data <- data_gen(n.t, N, p, rho, model.spec = 'GM', strata_num = strata.num)

# NEW CODE
setwd('~/Desktop/Stratified-SSL/R')
source('DataGeneration.R')
source('HelperFunctions.R')
set.seed(92047)
n_lab <- 400
n_unlab <- 20000
p <- 10
rho <- 0.4
num_strata <- 2
new_data <- DataGeneration(n_lab, n_unlab, p, rho, signal = c(1, 1, 0.5, 0.5),
                           model_specification = 'gaussian_mixture',
                           num_strata = num_strata)

# RUN TEST
test_result <- run_data_gen_test(old_data, new_data)
# Check if all tests passed
sum(unlist(test_result)) == length(test_result)
