setwd('~/Desktop/Stratified-SSL/R')
set.seed(92047)

## Step 1: Generate the data.
source('DataGeneration.R')

n_lab <- 400
n_unlab <- 20000
p <- 10
rho <- 0.4
model_specification <- 'outcome_correct_imputation_correct'
num_strata <- 2

my_data <- DataGeneration(n_lab, n_unlab, p, rho,
                         signal = c(1, 1, 0.5, 0.,5),
                         model_specification = model_specification,
                         num_strata = num_strata)

## Step 2: Get apparent parameter estimates.
