# Generate an example data (using the original code):
source('function.R')

N = 20000
p = 10
rho = 0.4
K.fold = 6
mod.spec = 'CC'
strata.num = 2
n.t = 200

data <- data_gen(n.t, N, p, rho, 
                 model.spec = mod.spec, strata_num = strata.num)

# Labled Data + Weights
Xt <- data$Xt
Yt <- data$Yt
Xv <- data$Xv
samp.prob <- data$samp.prob
Strata_ID <- data$S
St <- data$St
dat.all = rbind(Xt, Xv)

# Basis expansion
basis.x = ns.basis(dat.all, Strata_ID, nk = 3, basis.type = 'spline')



################# Modication of the DR approach #################
# functions in "DensityRatioRegression.R" and "DensityRatio"

source('DensityRatioRegression.R')
source('HelperFunctions.R')

basis_labeled <- basis.x[1:n.t,]
basis_unlabeled <- basis.x[-c(1:n.t),]
X_labeled <- Xt
X_unlabeled <- Xv
y <- Yt
samp_prob <- samp.prob

# Estimator for beta:

DR.fit <- DensityRatioRegression(basis_labeled, basis_unlabeled, X_labeled,
                                 X_unlabeled, y, samp_prob, lambda = NULL)
DR.fit$beta_dr

# Estimator for OMR and BS:

source('DensityRatioAccuracy.R')
DR.accuracy <- DensityRatioAccuracy(y, beta.dr = DR.fit$beta_dr, X_labeled, samp_prob, 
                                    DR.est = DR.fit$DR_est, c = 0.5)
# BS: 
DR.accuracy$mse.dr

# OMR:
DR.accuracy$ae.dr




################# Intrinsic efficient estimator #################

# See my newly added scripts "IntrinsicEfficientEst.R" for its implementation and
# scripts "IntrinsicSimulation.R" for generate the data for setup in Section S4.

source('HelperFunctions.R')
source('IntrinsicSimulation.R')
source('IntrinsicEfficientEst.R')


### Generate data 

n.t <- 400
rho <- 0.4
N <- 20000
p <- 3

# Setup 1

data <- IntrinsicData_set1(n.t, N, p, rho)
Xt <- data$Xt
Yt <- data$Yt
Xv <- data$Xv
samp.prob <- data$samp.prob
Strata_ID <- data$S
St <- data$St

# Basis expansion
dat.all = rbind(Xt, Xv)
basis.x = intri.basis(dat.all)


if (F){
  # Setup 2
  data <- IntrinsicData_set2(n.t, N, p, rho)
  Xt <- data$Xt
  Yt <- data$Yt
  Xv <- data$Xv
  samp.prob <- data$samp.prob
  Strata_ID <- data$S
  St <- data$St
  
  # Basis expansion
  dat.all = rbind(Xt, Xv)
  basis.x = intri.basis(dat.all)
}


# Solve for the plain SSL estimator as an initial estimator:

source('function')
glm.fit.ob = glm.fit.SS(basis.x, Xt, Xv, Yt, samp.prob, lambda0 = 0)
beta.ssl <- as.vector(glm.fit.ob$beta.ssl)

# Intrinsic estimate for beta:

# Index of the basis to be matched (i.e. X) in basis.x:
indx_mom <- c(1:2)
fit.result <- IntrinsicEfficientEstBeta(basis.x, Xt, Xv, Yt, samp.prob, 
                                        indx_mom, lambda0 = 0, theta_prelim = beta.ssl)
fit.result$theta


# Intrinsic estimate for Brier score:

est.bs <- IntrinsicEfficientEstBS(basis.x, Xt, Xv, Yt, samp.prob, indx_mom,
                                  lambda0 = 0, theta_prelim = beta.ssl)
est.bs$value


# Intrinsic estimate for OMR:

est.omr <- IntrinsicEfficientEstOMR(basis.x, Xt, Xv, Yt, samp.prob, indx_mom,
                                    lambda0 = 0, theta_prelim = beta.ssl, c = 0.5, h = NULL)
est.omr$value


