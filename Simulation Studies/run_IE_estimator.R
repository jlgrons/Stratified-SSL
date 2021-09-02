################################################################################
# We also provide an example for implementing the intrinsic efficiency estimator.
# This simulation study is provided in our Supplementary materials.

### Generate data 
set.seed(92047)
library('stratifiedSSL')
source('IE_helper_functions.R')
n.t <- 400
rho <- 0.4
N <- 20000
p <- 3

# Setup 1

data <- IntrinsicData_set1(n.t, N, p, rho)
X_labeled <- data$Xt
y <- data$Yt
X_unlabeled <- data$Xv
samp_prob <- data$samp.prob
Strata_ID <- data$S
S_labeled <- data$St

# Basis expansion
dat.all = rbind(X_labeled, X_unlabeled)
basis.x = intri.basis(dat.all)
basis_labeled = basis.x[1:n.t, ]
basis_unlabeled = basis.x[-(1:n.t), ]

if (F){
  # Setup 2
  data <- IntrinsicData_set2(n.t, N, p, rho)
  X_labeled <- data$Xt
  y <- data$Yt
  X_unlabeled <- data$Xv
  samp_prob <- data$samp.prob
  Strata_ID <- data$S
  S_labeled <- data$St
  
  # Basis expansion
  dat.all = rbind(X_labeled, X_unlabeled)
  basis.x = intri.basis(dat.all)
  basis_labeled = basis.x[1:n.t, ]
  basis_unlabeled = basis.x[-(1:n.t), ]
}


# Initial estimator.
beta_ssl <- SemiSupervisedRegression(basis_labeled, basis_unlabeled,
                                     X_labeled, X_unlabeled, y,
                                     samp_prob, lambda = 0)$beta_SSL

# Index of the basis to be matched (i.e. X) in basis.
indx_mom <- c(1:2)
IE_result <- IntrinsicEffEstBeta(basis_labeled, basis_unlabeled,
                                 X_labeled, X_unlabeled, y,
                                 samp_prob, indx_mom,
                                 lambda0 = 0, beta_ssl)
IE_result$theta


# Intrinsic estimate for Brier score:

IE_mse <- IntrinsicEfficientEstBS(basis_labeled, basis_unlabeled,
                                  X_labeled, X_unlabeled, y,
                                  samp_prob, indx_mom,
                                  lambda0 = 0, beta_ssl)
IE_mse$value


# Intrinsic estimate for OMR:

IE_omr <- IntrinsicEfficientEstOMR(basis_labeled, basis_unlabeled,
                                  X_labeled, X_unlabeled, y,
                                  samp_prob, indx_mom,
                                  lambda0 = 0, beta_ssl,
                                  threshold = 0.5, h = NULL)
IE_omr$value
