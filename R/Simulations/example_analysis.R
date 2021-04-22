## Step 1: Generate the data
n_lab <- 100
n_unlab <- 20000
p <- 10
rho <- 0.3
model_specification <- 'outcome_correct_imputation_correct'
num_strata <- 2

my_dat <- DataGeneration(n_lab, n_unlab, p, rho,
                         signal = c(1, 1, 0.5, 0.,5),
                         model_specification = model_specification,
                         num_strata = num_strata)
