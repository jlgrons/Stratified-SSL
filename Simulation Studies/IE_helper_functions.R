# Generate data (and the imputation basis) for Simulation in Section S4

# Setup (a) in Section S4:

IntrinsicData_set1 <- function(n.t, n.v, p, rho){
  
  N = n.t + n.v
  p = 3
  strata_num = 2
  Sigma0 = autocorr.mat(p = 3, rho = 0.4) 
  X0 = X.fun(N, Sigma = Sigma0);
  X_S = X0[,3]
  X0 = X0[,c(1:2)]
  S <- ifelse(I(X_S > 1), 1, 0)
  Y0 = S * (X0 %*% c(2, -2) + (-5) * X0[,1] * X0[,2]) + (1 - S) * (X0 %*% c(2, -2) + 5 * X0[,1] * X0[,2]) + rlogis(N)
  Y = I(Y0 > 0)
  
  n.each <- n.t / strata_num
  
  # Stratum 1 
  stratum.1 = which(S == 1);
  ind.1 = sample(stratum.1, size = n.each);
  
  # Stratum 2 
  stratum.2 = which(S == 0);
  ind.2 = sample(stratum.2, size = n.each);
  
  ### Labeled L & Unlabeled U datasets ###
  
  # Indices of labeled and unlabled data
  ind.t = c(ind.1, ind.2);  
  ind.v = setdiff(1:N, ind.t);  
  
  # Full data
  dat.N = cbind(Y, X0);
  
  # Labeled and Unlabed data sets
  
  dat.t = dat.N[ind.t, ]; 
  S.t = S[ind.t]
  dat.v = dat.N[ind.v, ]; 
  S.v = S[ind.v]
  S.all = c(S.t, S.v)
  
  # Relevant labeled data 
  Yt = dat.t[,1]; 
  Xt = dat.t[,-1]; 
  
  # Relevant unlabeled data
  Xv = dat.v[,-1]; 
  
  # Sampling probabilities
  ind.S1.v = (which(S > 0));
  ind.S2.v  = (which(S <= 0));
  N.1 = length(ind.S1.v)
  N.2 = length(ind.S2.v) 
  samp.prob = c(rep(n.each / N.1, n.each), rep(n.each / N.2, n.each));
  
  return(list(X0 = X0, Y = Y, S = S.all, St = S.t, Sv = S.v, Xt = Xt, 
              Xv = Xv, Yt = Yt, samp.prob = samp.prob))
  
}


# Setup (b) in Section S4:

IntrinsicData_set2 <- function(n.t, n.v, p, rho){
  
  N = n.t + n.v
  p = 3
  strata_num = 2
  Sigma0 = autocorr.mat(p = 3, rho = 0.4) 
  X0 = X.fun(N, Sigma = Sigma0);
  X_S = X0[,3]
  X0 = X0[,c(1:2)]
  S <- ifelse(I(X_S > 1), 1, 0)
  Y0 = X0 %*% c(1, -1) * (S + 1) + 1.5 * X0[,1] * X0[,2] + rlogis(N)
  Y = I(Y0 > 0)
  
  n.each <- n.t / strata_num
  
  # Stratum 1 
  stratum.1 = which(S == 1);
  ind.1 = sample(stratum.1, size = n.each);
  
  # Stratum 2 
  stratum.2 = which(S == 0);
  ind.2 = sample(stratum.2, size = n.each);
  
  ### Labeled L & Unlabeled U datasets ###
  
  # Indices of labeled and unlabled data
  ind.t = c(ind.1, ind.2);  
  ind.v = setdiff(1:N, ind.t);  
  
  # Full data
  dat.N = cbind(Y, X0);
  
  # Labeled and Unlabed data sets
  
  dat.t = dat.N[ind.t, ]; 
  S.t = S[ind.t]
  dat.v = dat.N[ind.v, ]; 
  S.v = S[ind.v]
  S.all = c(S.t, S.v)
  
  # Relevant labeled data 
  Yt = dat.t[,1]; 
  Xt = dat.t[,-1]; 
  
  # Relevant unlabeled data
  Xv = dat.v[,-1]; 
  
  # Sampling probabilities
  ind.S1.v = (which(S > 0));
  ind.S2.v  = (which(S <= 0));
  N.1 = length(ind.S1.v)
  N.2 = length(ind.S2.v) 
  samp.prob = c(rep(n.each / N.1, n.each), rep(n.each / N.2, n.each));
  
  return(list(X0 = X0, Y = Y, S = S.all, St = S.t, Sv = S.v, Xt = Xt, 
              Xv = Xv, Yt = Yt, samp.prob = samp.prob))
}



# Calculate the imputation basis:

intri.basis = function(X){
  
  basis <- cbind(X, X[,1] * X[,2])
  return(basis)
}



# Helper functions:

# Covariate generation
X.fun = function(NN, mu = rep(0,p), Sigma = Sigma0){mvrnorm(NN, mu, Sigma)};

# Autoregressive correlation structure
autocorr.mat <- function(p = 100, rho = 0.9) {
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}  




