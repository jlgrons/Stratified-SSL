num_sims <- 500

setwd('~/Desktop/Stratified-SSL/R/Original Code')
# Run with 200.

# Sample size of label data
n.t <- 200

# Matrix correlation
rho <- 0.4

# Sample size of unlabel data
N <- 20000
p <- 10
basis.type <- 'interact'
# Number of stratum
strata.num <- 4

source('response_main_IC.R')


# Run with 400.
setwd('~/Desktop/Stratified-SSL/R/Original Code')
# Sample size of label data
n.t <- 400

# Matrix correlation
rho <- 0.4

# Sample size of unlabel data
N <- 20000
p <- 10
basis.type <- 'interact'
# Number of stratum
strata.num <- 4

source('response_main_IC.R')
