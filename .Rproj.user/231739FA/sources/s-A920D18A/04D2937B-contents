#####  -----------------------------------------------------------------  #####
#####  Code for 'Efficient Estimation and Evaluation of Prediction Rules  #####
#####      in Semi-Supervised Settings under Stratified Sampling'         #####
#####  -----------------------------------------------------------------  #####


#########################################
### Necessary Functions and Libraries ###
#########################################


# Necessary packages
library(matrixStats)
library(stepPlr)
library(evd)
library(methods)
library(MASS)
library(glmnet)
library(randomForest)


# Useful functions for logistic model
dg.logit = function(xx){
  ddd <- exp(xx)/(exp(xx)+1)^2
  ddd[which(is.na(ddd))] = 0
  return(ddd)
}
logit = function(xx){log((xx)/(1-xx))};
g.logit = function(xx){1/(exp(- xx) + 1)}

# Covariate generation
X.fun = function(NN, mu = rep(0,p), Sigma = Sigma0){mvrnorm(NN, mu, Sigma)};

# Autoregressive correlation structure
autocorr.mat <- function(p = 100, rho = 0.9) {
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}

# Vector to matrix function
VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

# Computes the truncated cubic
trunc.cub = function(X, x){

  ## X: variable of interest
  ## x: knot location

  ((X>x)*(X-x))^3
}





########################
### Data Generation ###
#######################

data_gen <- function(n.t, n.v, p, rho, model.spec = 'CC', strata_num = 2){

  ## n.t: labeled data size
  ## n.v: unlabeled data size
  ## p: covariate dimension
  ## rho: correlation for covariates

  ## model.spec: model specification
  ## (both correct in Section 7 = 'CC',
  ## outcome model wrong, imputation model correct in Section 7 = 'IC',
  ## both wrong in Section 7 = 'II',
  ## outcome model wrong, imputation model correct in Section S4 = 'IC1',
  ## both wrong in Section S4 = 'II1',
  ## Gaussian mixture (GM) setting in Section S5 = 'GM')

  ## strata_num: number of stratum


  # Total data size
  N = n.t + n.v

  # Regression parameter
  b0 = c(1, 1, 0.5, 0.5); b0 = c(0, b0, rep(0, p - length(b0)));

  # Covariance for covariates
  Sigma0 = 3*autocorr.mat(p = p, rho = rho)

  ### Generate Covariates X0, Outcome Y, and Stratification variable S ###

  # Covariates
  X0 = X.fun(N,Sigma = Sigma0);

  # Linear predictor
  lin.pred = c(cbind(1,X0) %*% b0);

  # Outcome - Correct 'C' or Mis-specified Model otherwise
  if (model.spec == 'CC'){
    B = 0.5
    C = 0

    S1 <- I(X0[,1] + rnorm(N, 0, 1) < B)
    S3 <- I(X0[,3] + rnorm(N, 0, 1) < B)

    if (strata_num == 2){
      S <- ifelse(S1, 1, 0)
    }

    if (strata_num == 4){
      S <- rep(0, N)
      S[intersect(which(S1 == 1), which(S3 == 1))] <- 1
      S[intersect(which(S1 == 1), which(S3 == 0))] <- 2
      S[intersect(which(S1 == 0), which(S3 == 1))] <- 3
    }
    Y0 = lin.pred + rlogis(N)
    Y = I(Y0 > C)
  }

  if (model.spec == 'IC'){

    B = 0.5
    C = 0

    S1 <- I(X0[,1] + rnorm(N, 0, 1) < B)
    S3 <- I(X0[,3] + rnorm(N, 0, 1) < B)

    if (strata_num == 2){
      S <- ifelse(S1, 1, 0)
    }

    if (strata_num == 4){
      S <- rep(0, N)
      S[intersect(which(S1 == 1), which(S3 == 1))] <- 1
      S[intersect(which(S1 == 1), which(S3 == 0))] <- 2
      S[intersect(which(S1 == 0), which(S3 == 1))] <- 3
    }

    if (strata_num == 2){
      gamma.coef <- c(b0, c(0.5, 0, 0, 0.5), rep(0, 8), - 0.5, rep(0, 4), - 0.5)
      basis <- ns.basis(X0, S, 3, basis.type = 'interact')
    }
    if (strata_num == 4){
      gamma.coef <- c(b0, c(0.5, 0, 0, 0.5), rep(0, 8), - 0.5, rep(0, 4), - 0.5, 0, 0)
      basis <- ns.basis(X0, S, 3, basis.type = 'interact')
    }

    Y0 = cbind(1, basis) %*% gamma.coef + rlogis(N);
    Y = I(Y0 > C)
  }


  if (model.spec == 'II'){
    B = 0.5
    C = 0

    S1 <- I(X0[,1] + rnorm(N, 0, 1) < B)
    S3 <- I(X0[,3] + rnorm(N, 0, 1) < B)

    if (strata_num == 2){
      S <- ifelse(S1, 1, 0)
    }

    if (strata_num == 4){
      S <- rep(0, N)
      S[intersect(which(S1 == 1), which(S3 == 1))] <- 1
      S[intersect(which(S1 == 1), which(S3 == 0))] <- 2
      S[intersect(which(S1 == 0), which(S3 == 1))] <- 3
    }
    c0 = c(0, -2, rep(0,3), -3, -3, 0, 0, rep(0,2))
    Y0 = lin.pred + X0[,1]^2 + X0[,3]^2 + rgumbel(N, -2, 0.3)*exp(cbind(1,X0)%*% c0);
    Y = I(Y0 > C)
  }

  if (model.spec == 'IC1'){
    B = 1.5
    C = 5

    S1 <- I(X0[,1] + X0[,2] + rnorm(N, 0, 1) > B)
    if (strata_num == 2){
      S <- ifelse(S1, 1, 0)
    }

    gamma.coef <- c(b0, c(0.5, 0, 0, 0.5), rep(0, 8), - 0.5, rep(0, 4), - 0.5)
    basis <- ns.basis(X0, S, 3, basis.type = 'interact')

    Y0 = cbind(1, basis) %*% gamma.coef * S +
      (0.8 * cbind(1, basis) %*% gamma.coef - C) * (1 - S) + rlogis(N);
    Y = I(Y0 > 1)
  }

  if (model.spec == 'II1'){
    B = 1.5
    C = 5

    S1 <- I(X0[,1] + X0[,2] + rnorm(N, 0, 1) > B)
    if (strata_num == 2){
      S <- ifelse(S1, 1, 0)
    }

    c0 = c(0, -2, rep(0,3), -3, -3, 0, 0, rep(0,2))
    Y0 = (lin.pred + X0[,1]^2 + X0[,3]^2) * S + (0.8 * lin.pred - C)* (1 - S) +
      rgumbel(N, -2, 0.3)*exp(cbind(1,X0)%*% c0);
    Y = I(Y0 > 1)
  }


  if (model.spec == 'GM'){
    B = 0.5
    mu_diff <- c(0.2, -0.2, 0.2, -0.2, 0.2, -0.2, 0.1, -0.1, rep(0, p - 8))
    Sigma_diff <- autocorr.mat(p, rho = 0.3) - diag(rep(0.4, p)) + matrix(0.2, p, p)

    #S1 <- I(X0[,1] + rnorm(N, 0, 1) < B)
    #S3 <- I(X0[,3] + rnorm(N, 0, 1) < B)

    #if (strata_num == 2){
    #  S <- ifelse(S1, 1, 0)
   # }

    #if (strata_num == 4){
    #  S <- rep(0, N)
    #  S[intersect(which(S1 == 1), which(S3 == 1))] <- 1
    #  S[intersect(which(S1 == 1), which(S3 == 0))] <- 2
    #  S[intersect(which(S1 == 0), which(S3 == 1))] <- 3
    #}

    Y = rbinom(N, 1, 0.5)
    X0 = matrix(0, N, p)

    mu1 = rep(0, p)
    Sigma1 = autocorr.mat(p = p, rho = rho)
    X0[which(Y == 0), ] = X.fun(length(which(Y == 0)), mu1, Sigma1)

    mu2 = mu1 + mu_diff
    Sigma2 = Sigma1 + Sigma_diff
    X0[which(Y == 1), ] = X.fun(length(which(Y == 1)), mu2, Sigma2)

    S1 <- I(X0[,1] + rnorm(N, 0, 1) < B)
    S3 <- I(X0[,3] + X0[,4] + rnorm(N, 0, 1) > B)

    if (strata_num == 2){
      S <- ifelse(S1, 1, 0)
    }

    if (strata_num == 4){
      S <- rep(0, N)
      S[intersect(which(S1 == 1), which(S3 == 1))] <- 1
      S[intersect(which(S1 == 1), which(S3 == 0))] <- 2
      S[intersect(which(S1 == 0), which(S3 == 1))] <- 3
    }
    X0[which(Y == 1), 3] <- X0[which(Y == 1), 3] + 0.12 * X0[which(Y == 1), 3]^3
    X0[which(Y == 1), 4] <- X0[which(Y == 1), 4] + 0.12 * X0[which(Y == 1), 4]^3

    X0[which(Y == 1), 7] <- X0[which(Y == 1), 7] + 0.12 * X0[which(Y == 1), 7]^3
    X0[which(Y == 1), 8] <- X0[which(Y == 1), 8] + 0.12 * X0[which(Y == 1), 8]^3

  }

  # Size of strata

  n.each <- n.t / strata_num

  ##ind_lst <- vector('list', strata_num)

  if (strata_num == 2){

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

    ## Stratified sampling

    dat.t = dat.N[ind.t, ];
    S.t = S[ind.t]
    dat.v = dat.N[ind.v, ];
    S.v = S[ind.v]
    S.all = c(S.t, S.v)

    ## Random sampling (only used under our settings in Section S4)

    dat.r = dat.N[1:n.t, ]
    Yr = dat.r[,1];
    Xr = dat.r[,-1];
    Xvr = dat.N[-c(1:n.t), -1];

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
  }

  if (strata_num == 4){

    # Stratum 1
    stratum.1 = which(S == 0);
    ind.1 = sample(stratum.1, size = n.each);

    # Stratum 2
    stratum.2 = which(S == 1);
    ind.2 = sample(stratum.2, size = n.each);

    # Stratum 3
    stratum.3 = which(S == 2);
    ind.3 = sample(stratum.3, size = n.each);

    # Stratum 4
    stratum.4 = which(S == 3);
    ind.4 = sample(stratum.4, size = n.each);

    ### Labeled L & Unlabeled U datasets ###

    # Indices of labeled and unlabled data
    ind.t = c(ind.1, ind.2, ind.3, ind.4);
    ind.v = setdiff(1:N, ind.t);

    # Full data
    dat.N = cbind(Y, X0);

    # Labeled and Unlabed data sets
    dat.t = dat.N[ind.t, ];
    S.t = S[ind.t]
    dat.v = dat.N[ind.v, ];
    S.v = S[ind.v]
    S.all = c(S.t, S.v)

    ## Random sampling (only used under our settings in Section S4)

    dat.r = dat.N[1:n.t, ]
    Yr = dat.r[,1];
    Xr = dat.r[,-1];
    Xvr = dat.N[-c(1:n.t), -1];

    # Relevant labeled data
    Yt = dat.t[,1];
    Xt = dat.t[,-1];

    # Relevant unlabeled data
    Xv = dat.v[,-1];

    # Sampling probabilities
    ind.S1.v = (which(S == 0));
    ind.S2.v  = (which(S == 1));
    ind.S3.v = (which(S == 2));
    ind.S4.v  = (which(S == 3));
    N.1 = length(ind.S1.v)
    N.2 = length(ind.S2.v)
    N.3 = length(ind.S3.v)
    N.4 = length(ind.S4.v)
    samp.prob = c(rep(n.each/N.1, n.each), rep(n.each/N.2, n.each),
                  rep(n.each/N.3, n.each), rep(n.each/N.4, n.each));

  }

  return(list(X0 = X0, Y = Y, S = S.all, St = S.t, Sv = S.v, Xt = Xt,
              Xv = Xv, Yt = Yt, samp.prob = samp.prob, b0 = b0,
              Xr = Xr, Yr = Yr, Sr = S, Xvr = Xvr))

}


##############
### Method ###
##############

# Computes the natural spline basis
ns.basis = function(X, S, nk, basis.type = NULL){

  ## X: covariate matrix
  ## S: the stratum ID
  ## nk: number of knots
  ## basis.type: type of basis (choices under different settings)

  if (basis.type == 'spline'){
    X = as.matrix(X); basis = NULL;

    for(i in 1:ncol(X)){
      X_i = X[,i];

      # quantiles to determine appropriate knots
      knots <- quantile(X_i, seq(0, 1, length = nk));

      # changes quantiles if there aren't enough unique values
      j = 0;
      while(length(unique(knots)) != nk){
        j = j + 1;
        knots <- unique(quantile(X_i, seq(0, 1, length = nk  + j)))
      }

      # compute the natural spline basis
      d_k <- (trunc.cub(X_i, knots[nk-1])) /(knots[nk] - knots[nk-1]);
      evals <- sapply(1:(nk-2), function(ii){
        d_i <- (trunc.cub(X_i, knots[ii])) /(knots[nk] - knots[ii]);
        basis.new <- d_i - d_k})

      # bind together and include original variable
      basis = cbind(basis, cbind(X_i, evals))

    }
    basis.S <- c()
    for (k in 1:(length(unique(S)) - 1)) {
      basis.S <- cbind(basis.S, ifelse(S == (k - 1), 1, 0))
    }
    basis <- cbind(basis, basis.S)
  }

  if (basis.type == 'interact'){
    basis <- X
    p <- length(X[1,])
    for (j in 2:p){
      basis <- cbind(basis, X[,1] * X[,j])
    }
    for (j in 3:p){
      basis <- cbind(basis, X[,2] * X[,j])
    }
    basis.S <- c()
    for (k in 1:(length(unique(S)) - 1)) {
      basis.S <- cbind(basis.S, ifelse(S == (k - 1), 1, 0))
    }
    basis <- cbind(basis, basis.S)
  }

  if (basis.type == 'IC1'){
    basis <- X
    p <- length(X[1,])
    for (j in 2:p){
      basis <- cbind(basis, X[,1] * X[,j])
    }
    for (j in 3:p){
      basis <- cbind(basis, X[,2] * X[,j])
    }
    basis.S <- S
    for (k in 1:length(basis[1, ])) {
      basis.S <- cbind(basis.S, S * basis[,k])
    }
    basis <- cbind(basis, basis.S)
  }

  if (basis.type == 'II1'){
    X = as.matrix(X); basis = NULL;

    for(i in 1:ncol(X)){
      X_i = X[,i];

      # quantiles to determine appropriate knots
      knots <- quantile(X_i, seq(0, 1, length = nk));

      # changes quantiles if there aren't enough unique values
      j = 0;
      while(length(unique(knots)) != nk){
        j = j + 1;
        knots <- unique(quantile(X_i, seq(0, 1, length = nk  + j)))
      }

      # compute the natural spline basis
      d_k <- (trunc.cub(X_i, knots[nk-1])) /(knots[nk] - knots[nk-1]);
      evals <- sapply(1:(nk-2), function(ii){
        d_i <- (trunc.cub(X_i, knots[ii])) /(knots[nk] - knots[ii]);
        basis.new <- d_i - d_k})

      # bind together and include original variable
      basis = cbind(basis, cbind(X_i, evals))

    }
    basis <- cbind(basis, S)
    for (k in 1:length(X[1,])) {
      basis <- cbind(basis, X[,k] * S)
    }

  }

  return(basis)
}



##############################
### SSL GLMs via Basis Exp ###
##############################

# fits the penalized ridge model for semi-supervised estimate
my.ridge = function(x, y, weights = NULL, lambda0 = 1e-04){

  ## x: predictors
  ## y: outcome
  ## weights: sampling or resampling weights

  if(is.null(weights)){weights = rep(1, length(y))}

  # multiple lambdas to make sure convergence happens
  gamma = tryCatch((coef(glmnet(x, y, weights = weights, alpha = 0,
                                lambda = seq(lambda0, 100*lambda0, length.out = 100),
                                family = 'binomial'))),
                   error = function(e) rep(NA, ncol(x) + 1 ))

  # take the value at lambda0
  gamma = gamma[, ncol(gamma)]
  return(gamma)
}


# fits the glm using unlabeled and labeled data & computes corresponding supervised estimator
glm.fit.SS = function(basis.x, Xt, Xv, Yt, samp.prob, lambda0 = NULL){

  ## basis.x: basis matrix with first n.t rows corr to labeled data and unlabeled rows n.t+1 to n.t+n.v
  ## Xt: covariates for labeled
  ## Xt: covariates for labeled
  ## Yt: labels
  ## samp.prob: vector of sampling weights
  ## lambda0: ridge parameter for SS estimator

  n.t = length(Yt); pp = ncol(basis.x); p = ncol(Xt); ind.lab = 1:n.t;

  if(is.null(lambda0)){lambda0 = log(pp)/n.t^1.5};

  weights = 1/samp.prob/mean(1/samp.prob)

  # step 1: basis function regression
  gamma = my.ridge(basis.x[ind.lab, ], Yt, weights = weights, lambda = lambda0);

  # step 2: fit the SS GLM
  imps.basis = g.logit(cbind(1, basis.x) %*% gamma);
  dat.all = rbind(Xt, Xv);
  beta.ssl = tryCatch(glm(imps.basis ~ dat.all,
                          family = 'binomial')$coeff, error = function(e) rep(NA, p+1))

  # fit supervised estimate for comparison
  beta.sl = tryCatch(glm(Yt ~ Xt, family = 'binomial',
                         weights = weights)$coeff, error = function(e) rep(NA, p+1));

  beta.naive = tryCatch(glm(Yt ~ Xt, family = 'binomial')$coeff,
                        error = function(e) rep(NA, p+1));

  # Fit the density ratio estimator

  n.v <- length(basis.x[,1])
  basis.phi <- cbind(rep(1, n.v), basis.x)
  phi.t <- basis.phi[1:n.t, ]
  phi.v <- basis.phi[- c(1:n.t), ]
  dim.basis <- length(basis.phi[1, ])

  phiT.phi <- t(basis.phi) %*% basis.phi / (n.t + n.v)
  E.phi <- t(phi.t) %*% weights / n.t
  theta.ratio <- solve(phiT.phi + diag(rep(0.01, ncol(phiT.phi)))) %*% (rowMeans(t(phi.v)) - E.phi)
  density.ratio <- exp(phi.t %*% theta.ratio)
  weights.dr <- density.ratio * weights

  beta.dr <- tryCatch(glm(Yt~Xt, family = 'binomial',
                          weights = weights.dr)$coeff, error = function(e) rep(NA, p+1));

  u.dr <- diag(as.vector(Yt - g.logit(cbind(1, Xt) %*% beta.dr))) %*% cbind(1, Xt)
  E.uT.phi <- t(u.dr) %*% diag(as.vector(weights)) %*% phi.t / sum(weights)
  proj.coef.dr <- E.uT.phi %*% solve(phiT.phi + diag(rep(0.01, ncol(phiT.phi))))
  proj.dr <- diag(as.vector(weights)) %*% phi.t %*% t(proj.coef.dr)

  return(list(beta.ssl = beta.ssl, beta.sl = beta.sl, gamma = gamma,
              beta.dr = beta.dr, proj.dr = proj.dr,
              beta.naive = beta.naive))
}


#################################
### Cross-Validated Residuals ###
#################################

# Computes cross-validated (CV) residuals based on supervised and SS ests of the regression parameter
# as well as the regression parameter from the imputation model
resids.cv = function(basis.x, Xt, Xv, Yt, samp.prob, K.fold, lambda0 = NULL){

  ## basis.x: basis matrix with first n.t rows corr to labeled data and unlabeled rows n.t+1 to n.t+n.v
  ## Xt: covariates for labeled
  ## Xv: covariates for unlabeled
  ## Yt: labels
  ## samp.prob: vector of sampling weights
  ## K.fold: # of folds for CV
  ## lambda0: ridge parameter for SS estimator

  n.t = nrow(Xt); ind.cv = split(1:n.t, sample(rep(1:K.fold, floor(n.t/K.fold))));
  qq = length(ind.cv); ind.lab = 1:n.t; pp = ncol(basis.x);

  if(is.null(lambda0)){lambda0 = log(pp)/(floor((K.fold-1)*n.t/K.fold))^1.5};

  # K estimates of supervised and SS beta & gamma with kth fold removed
  cv.ests = lapply(1:qq, function(kk) {inds.t = as.vector(unlist(ind.cv[-kk]));
  Yt.t = Yt[inds.t]; samp.prob.t = samp.prob[inds.t];
  basis.cv = rbind(basis.x[inds.t,], basis.x[-ind.lab,]);
  glm.fit.SS(basis.cv,  Xt[inds.t,], Xv, Yt.t, samp.prob.t, lambda0 = lambda0)});

  beta.ssl.cv = sapply(cv.ests, "[[", 1);
  beta.sl.cv = sapply(cv.ests, "[[", 2);
  gamma.cv = sapply(cv.ests, "[[", 3);
  beta.dr.cv = sapply(cv.ests, "[[", 4);

  # Basis corresponding to the labeled set
  basis.lab = basis.x[ind.lab, ];

  # K sets of residuals based on the kth fold and the beta with the kth fold removed
  res.cv = lapply(1:qq,function(kk){inds.v = as.vector(unlist(ind.cv[kk]));
  beta.ssl.cv.tmp = beta.ssl.cv[,kk]; beta.sl.cv.tmp = beta.sl.cv[,kk];
  gamma.cv.tmp = gamma.cv[,kk]; beta.dr.cv.tmp = beta.dr.cv[,kk];
  samp.prob.v = samp.prob[inds.v];
  Yt.v = Yt[inds.v]; Xt.v = Xt[inds.v,]; basis.v = basis.lab[inds.v, ];
  pred.b.ssl = g.logit(cbind(1, Xt.v) %*% beta.ssl.cv.tmp);
  pred.b.sl = g.logit(cbind(1, Xt.v) %*% beta.sl.cv.tmp);
  pred.gamma = g.logit(cbind(1, basis.v) %*% gamma.cv.tmp);
  pred.b.dr = g.logit(cbind(1, Xt.v) %*% beta.dr.cv.tmp);
  resids = list(beta.ssl = (Yt.v - pred.b.ssl)*(1/samp.prob.v),
                beta.sl = (Yt.v - pred.b.sl)*(1/samp.prob.v),
                gamma = (Yt.v - pred.gamma)*(1/samp.prob.v),
                beta.dr = (Yt.v - pred.b.dr)*(1/samp.prob.v))
  })

  inds.ord = order(unlist(ind.cv)); wgt = mean(1/samp.prob);

  # Reorder the residuals and divide by the mean of the weights
  resids.beta.ssl = c(unlist(sapply(res.cv, "[[", 1)))[inds.ord]/wgt;
  resids.beta.sl = c(unlist(sapply(res.cv, "[[", 2)))[inds.ord]/wgt;
  resids.gamma = c(unlist(sapply(res.cv, "[[", 3)))[inds.ord]/wgt;
  resids.beta.dr = c(unlist(sapply(res.cv, "[[", 4)))[inds.ord]/wgt;


  return(list(resids.beta.ssl = resids.beta.ssl, resids.beta.sl = resids.beta.sl,
              resids.gamma = resids.gamma, resids.beta.dr = resids.beta.dr,
              beta.ssl.cv = beta.ssl.cv, beta.sl.cv = beta.sl.cv,
              gamma.cv = gamma.cv, beta.dr.cv = beta.dr.cv))
}

# computes a standard error estimate using the influence function of
# solution to a simple estimating equation
se.est = function(bhat, Yt, Xt, Xv, resid){

  ## bhat: parameter estimate
  ## Xt: covariates for labeled
  ## Xv: covariates for unlabeled
  ## Yt: labels
  ## resids: vector of CV/non-CV residuals to compute score

  Xt.1 = cbind(1,Xt); Xv.1 = cbind(1, rbind(Xt,Xv));
  n.v = nrow(Xv.1); n.t = nrow(Xt.1);

  S =  t(Xt.1) %*% (Xt.1 * c(resid)^2)/n.t;
  ddd = c(dg.logit(Xv.1%*%bhat))
  ddd[which(is.na(ddd))] = 0
  I = t(Xv.1) %*% (Xv.1 * ddd)/n.v;

  A = solve(I);
  vars =  A %*% S %*% A;
  return(list(se = sqrt(diag(vars)/n.t), A = A))
}


se.est.dr <- function(bhat, Yt, Xt, Xv, resid, proj.dr){
  Xt.1 = cbind(1, Xt); Xv.1 = cbind(1, rbind(Xt, Xv));
  n.v = nrow(Xv.1); n.t = nrow(Xt.1);

  u.dr <- diag(as.vector(resid)) %*% Xt.1
  S =  t(u.dr - proj.dr) %*% (u.dr - proj.dr) /n.t;

  ddd = c(dg.logit(Xv.1 %*% bhat))
  ddd[which(is.na(ddd))] = 0

  I = t(Xv.1) %*% (Xv.1 * ddd)/n.v;

  A = solve(I);
  vars =  A %*% S %*% A;
  return(sqrt(diag(vars)/n.t))

}


#############################################
### Component-wise minimum var. estimator ###
#############################################

# computes the component-wise minimum variance estimator based on the CV resids
min.var.est = function(beta.ssl, beta.sl, resids.gamma, resids.beta.sl,
                       Xt, Xv, eps = NULL){
  Xv.1 = cbind(1, rbind(Xt,Xv)); Xt.1 = cbind(1,Xt);
  n.v = nrow(Xv.1); n.t = nrow(Xt.1); ones = c(1,1);
  p = ncol(Xt); weight = matrix(NA, p +1, 2);

  if(is.null(eps)){eps = rep(0, p+1)};

  # compute minimum var estimator (note the residuals have been divided by mean of weights)

  A = solve(t(Xv.1) %*% (Xv.1*c(dg.logit(Xv.1%*%beta.ssl))))*n.v;
  T_1 = A %*% t(Xt.1*c(resids.gamma));
  T_2 = A %*% t(Xt.1*c(resids.beta.sl));

  for (i in 1:(p + 1)){
    mat = cbind(T_1[i,], T_2[i,]);
    cov = solve(t(mat)%*%mat  + eps[i]*diag(2))/n.t;
    w2= t(ones) %*% cov %*% ones
    w1 = t(ones) %*% cov
    weight[i,] = w1/c(w2)
    if (NA %in% weight[i,]){
      weight[i,] <- c(1, 0)
    }

  }

  # min var estimator
  beta.ssl.w = beta.ssl*weight[,1] + weight[,2]*beta.sl

  # std error of min var
  A = solve(t(Xv.1) %*% (Xv.1*c(dg.logit(Xv.1%*%beta.ssl.w))))*n.v;
  w.beta.ssl = diag(weight[,1]); w.beta.sl =  diag(weight[,2]);
  resids.w = (w.beta.ssl %*% t(Xt.1*c(resids.gamma)) +
                w.beta.sl %*% t(Xt.1*resids.beta.sl))/n.t
  se.beta.w = sqrt(diag(A %*% (resids.w %*% t(resids.w)) %*% A))


  return(list(beta = beta.ssl.w, weight= weight, se.est = se.beta.w))
}


################################################################################
##### Cross-validated estimator of the model accuracy parameter (BS and OMR) ###
################################################################################

# computes the CV estimates for model evaluation parameters
model.eval.cv = function(basis.x, Xv, Xt, Yt, samp.prob,
                         W, K.fold = 5, rep = 10, c = 0.5, lambda0 = NULL){

  ## basis.x: basis exp matrix
  ## Xt: covariates for labeled
  ## Xv: covariates for unlabeled
  ## Yt: labels
  ## samp.prob: sampling probability
  ## W: weight for the minimum variance estimator
  ## K.fold: number of folds for CV
  ## rep: number of reps for CV
  ## c: cutoff for the overall misclassification rate
  ## lambda0: ridge parameter for the imputation model

  pp = ncol(basis.x); n.t = nrow(Xt); ind.lab = 1:n.t;
  data.all = rbind(Xt,Xv);

  mse.cv.ssl = matrix(NA, nrow = K.fold, ncol = rep);
  ae.cv.ssl = matrix(NA, nrow = K.fold, ncol = rep);

  mse.cv.sl = matrix(NA, nrow = K.fold, ncol = rep);
  ae.cv.sl = matrix(NA, nrow = K.fold, ncol = rep);

  mse.cv.naive = matrix(NA, nrow = K.fold, ncol = rep);
  ae.cv.naive = matrix(NA, nrow = K.fold, ncol = rep);

  mse.cv.dr = matrix(NA, nrow = K.fold, ncol = rep);
  ae.cv.dr = matrix(NA, nrow = K.fold, ncol = rep);


  if(is.null(lambda0)){lambda0 = log(pp)/(floor((K.fold-1)*n.t/K.fold))^1.5};

  for(j in 1:rep){
    set.seed(j)
    ind.cv = split(1:n.t, sample(rep(1:K.fold, floor(n.t/K.fold))));

    for(i in 1:K.fold){

      inds.v = as.vector(unlist(ind.cv[i]));
      inds.t = as.vector(setdiff(as.vector(unlist(ind.cv)), inds.v))

      wg.v = samp.prob[inds.v]; wg.t = samp.prob[inds.t]

      Yt.v = Yt[inds.v]; Xt.v = Xt[inds.v,]; basis.v = basis.x[inds.v, ];
      Yt.t = Yt[inds.t]; Xt.t = Xt[inds.t,]; basis.t = basis.x[inds.t, ]

      beta.tmp = glm.fit.SS(rbind(basis.t, basis.x[-ind.lab, ]),
                            Xt.t, Xv, Yt.t, wg.t, lambda0)

      beta.ssl.t =  beta.tmp$beta.ssl
      beta.sl.t = beta.tmp$beta.sl
      gamma.t =  beta.tmp$gamma
      beta.ssl.w.t = W*beta.ssl.t + (1-W)*beta.sl.t;
      beta.dr.t = beta.tmp$beta.dr

      lp.v = g.logit(cbind(1,Xt.v)%*%beta.ssl.w.t);
      lp.v.ind = ifelse(I(lp.v > c), 1, 0);

      imps.v = cbind(1, basis.v)%*%gamma.t;
      refit.pv.1 <- glm(Yt.v~cbind(lp.v), offset = imps.v,
                        weights = 1/wg.v/mean(1/wg.v), family = 'binomial')$coeff
      refit.pv.2 <-glm(Yt.v~cbind(lp.v.ind), offset = imps.v,
                       weights = 1/wg.v/mean(1/wg.v), family = 'binomial')$coeff

      lp.u = g.logit(cbind(1,data.all)%*%beta.ssl.w.t);
      lp.u.ind = I(lp.u > c);

      imps.pe.1 = g.logit(cbind(1,lp.u)%*% refit.pv.1 + cbind(1, basis.x)%*%gamma.t);
      imps.pe.2 = g.logit(cbind(1,lp.u.ind)%*% refit.pv.2 + cbind(1, basis.x)%*%gamma.t);

      mse.cv.ssl[i,j] = mean(imps.pe.1 + (lp.u - 2*imps.pe.1)*lp.u)
      ae.cv.ssl[i,j] = mean(imps.pe.2 + (lp.u.ind - 2*imps.pe.2)*lp.u.ind)

      lp.v.sl = g.logit(cbind(1,Xt.v) %*% beta.sl.t);
      lp.v.sl.ind = I(lp.v.sl > c);

      mse.cv.sl[i,j] = mean((Yt.v  - lp.v.sl)^2 * 1/wg.v)/mean(1/wg.v)
      ae.cv.sl[i,j] = mean(abs(Yt.v - lp.v.sl.ind) * 1/wg.v)/mean(1/wg.v)

      mse.cv.naive[i,j] = mean((Yt.v  - lp.v.sl)^2)
      ae.cv.naive[i,j] = mean(abs(Yt.v - lp.v.sl.ind))

      lp.v.dr = g.logit(cbind(1, Xt.v) %*% beta.dr.t);
      lp.v.dr.ind = I(lp.v.dr > c);

      mse.cv.dr[i,j] = mean((Yt.v  - lp.v.dr)^2 * 1 / wg.v) / mean(1 / wg.v)
      ae.cv.dr[i,j] = mean(abs(Yt.v - lp.v.dr.ind) * 1 / wg.v) / mean(1 / wg.v)

    }
  }

  return(list(mse.ssl = mean(mse.cv.ssl, na.rm = T),
              ae.ssl = mean(ae.cv.ssl, na.rm = T),
              mse.sl = mean(mse.cv.sl, na.rm = T),
              ae.sl = mean(ae.cv.sl, na.rm = T),
              mse.dr = mean(mse.cv.dr, na.rm = T),
              ae.dr = mean(ae.cv.dr, na.rm = T),
              mse.naive = mean(mse.cv.naive, na.rm = T),
              ae.naive = mean(ae.cv.naive, na.rm = T)))

}



###########################################################################
####### Apparent estimator of the model accuracy parameter (BS and OMR) ###
###########################################################################

# computes the apparent estimates for model evaluation parameters
model.eval.ap = function(Yt, beta.sl, beta.ssl, gamma, beta.dr,
                         Xt, Xv, basis.x, samp.prob, V = NULL, c = 0.5){

  ## beta.sl: supervised beta
  ## beta.ssl: semi-supervised beta
  ## gamma: estimated parameter for imputation model
  ## Xt: covariates for labeled
  ## Xv: covariates for unlabeled
  ## Yt: labels
  ## basis.x: basis exp matrix
  ## samp.prob: sampling probability
  ## V: optional weights for resampling
  ## c: cutoff for the overall misclassification rate

  if(is.null(V)){V = rep(1, length(Yt))};
  weight = V/samp.prob/mean(V/samp.prob);
  dat.all = rbind(Xt, Xv);

  lp.v.ssl = g.logit(cbind(1,dat.all)%*%beta.ssl);
  lp.v.ssl.ind = I(lp.v.ssl > c);

  lp.t.ssl = g.logit(cbind(1,Xt)%*%beta.sl);
  lp.t.ssl.ind = I(lp.t.ssl > c);

  lp.t.dr = g.logit(cbind(1,Xt)%*%beta.dr);

  lp.t.dr.ind = I(lp.t.dr > c);

  n.t = length(Yt); ind.l = (1:n.t);
  imps.t = cbind(1, basis.x[ind.l,])%*%gamma

  refit.p1 <- glm(Yt~cbind(lp.t.ssl), offset = imps.t, family = 'binomial', weights = weight)$coeff;
  refit.p2 <- glm(Yt~cbind(lp.t.ssl.ind), offset = imps.t, family = 'binomial', weights = weight)$coeff;
  imps.pe1 = g.logit(cbind(1,lp.v.ssl)%*%refit.p1 +  cbind(1, basis.x)%*%gamma);
  imps.pe2 = g.logit(cbind(1,lp.v.ssl.ind)%*%refit.p2 +  cbind(1, basis.x)%*%gamma);

  mse.ssl = mean(imps.pe1 + (lp.v.ssl - 2*imps.pe1)*lp.v.ssl)
  ae.ssl = mean(imps.pe2 + (lp.v.ssl.ind - 2*imps.pe2)*lp.v.ssl.ind)

  mse.sl = mean((Yt - lp.t.ssl)^2 * weight)
  ae.sl = mean(abs(Yt - lp.t.ssl.ind) * weight)

  mse.naive = mean((Yt - lp.t.ssl)^2)
  ae.naive = mean(abs(Yt - lp.t.ssl.ind))

  mse.dr = mean((Yt - lp.t.dr)^2 * weight)
  ae.dr = mean(abs(Yt - lp.t.dr.ind) * weight)

  return(list(mse.ssl = mse.ssl, ae.ssl = ae.ssl,
              mse.sl = mse.sl, ae.sl = ae.sl, mse.dr = mse.dr,
              ae.dr = ae.dr, mse.naive = mse.naive, ae.naive = ae.naive))
}

#######################################################
### SE Est for model accuracy parameters (BS & OMR) ###
#######################################################

# computes the standard error estimates for model evaluation parameters
model.eval.se = function(b, Yt, Xt, Xv, basis.x, samp.prob, W, beta.sl, beta.ssl.w,
                         beta.dr, resids.beta.sl, resids.gamma, resids.beta.dr,
                         A, c = 0.5){

  ## b: number of resamples
  ## Xt: covariates for labeled
  ## Xv: covariates for unlabeled
  ## Yt: labels
  ## basis.x: basis exp matrix
  ## samp.prob: sampling probability
  ## W: weight for the minimum variance estimator
  ## beta.sl: supervised beta
  ## beta.ssl.w: semi-supervised beta (minVar)
  ## resids.beta.sl: CV resids for supervised beta
  ## resids.gamma: CV resids for imputation parameter
  ## A: information matrix based on semi-supervisd beta (minVar)
  ## c: cutoff for the overall misclassification rate

  wgt = sapply(1:b, function(kk) 4*rbeta(n.t, 0.5, 1.5))

  resids.beta.sl.p = ((wgt-1)*resids.beta.sl)
  resids.gamma.p = ((wgt-1)*resids.gamma)
  resids.beta.dr.p = (wgt-1)*resids.beta.dr
  proj.dr.p = t(proj.dr) %*% (wgt-1)

  Xt.1 = cbind(1, Xt)
  T_1.p = A %*% t(Xt.1)%*%resids.gamma.p/n.t;
  T_2.p = A %*% t(Xt.1)%*%resids.beta.sl.p/n.t;

  T_1.dr.p = A %*% t(Xt.1) %*% resids.beta.dr.p / n.t
  T_2.dr.p = A %*% proj.dr.p / n.t

  beta.ssl.w.p = beta.ssl.w + diag(W[,1]) %*% T_1.p + diag(W[,2]) %*% T_2.p;
  beta.sl.p = beta.sl + T_2.p
  beta.dr.p = beta.dr + T_1.dr.p - T_2.dr.p

  ind.lab = 1:n.t;

  gamma.p = sapply(1:b, function(kk){
    my.ridge(basis.x[ind.lab, ], Yt, weights = wgt[,kk]/samp.prob/mean(wgt[,kk]/samp.prob),
             lambda = log(ncol(basis.x))/(n.t^1.5))
  })

  perts = lapply(1:b, function(kk){model.eval.ap(Yt, beta.sl.p[,kk], beta.ssl.w.p[,kk], gamma.p[,kk],
                                                 beta.dr.p[,kk], Xt, Xv, basis.x, samp.prob, V = wgt[,kk])})

  ssl.pert.mse = sapply(perts, "[[", 1);
  ssl.pert.ae = sapply(perts, "[[", 2);

  sl.pert.mse = sapply(perts, "[[", 3);
  sl.pert.ae = sapply(perts, "[[", 4);

  dr.pert.mse = sapply(perts, "[[", 5);
  dr.pert.ae = sapply(perts, "[[", 6);

  return(list(ssl.pert.mse = ssl.pert.mse, sl.pert.mse = sl.pert.mse,
              dr.pert.mse = dr.pert.mse, ssl.pert.ae = ssl.pert.ae,
              sl.pert.ae = sl.pert.ae, dr.pert.ae = dr.pert.ae,
              gamma.p = gamma.p, beta.ssl.w.p = beta.ssl.w.p,
              beta.sl.p = beta.sl.p, beta.dr.p = beta.dr.p))

}


######################################################################
############ Extract direction of the principle component ############
############ Used for the ML settings in Section S5  #################
######################################################################


extract.pc <- function(X, out.num = 5, basis.num = 10,
                       V.first = NULL, V.second = NULL){
  p <- length(X[1,])
  X_sec_mat <- c()
  for (j in 1:p){
    for (k in (j:p)){
      X_sec_mat <- cbind(X_sec_mat, X[,j] * X[,k])
    }
  }

  if (is.null(V.first)){
    svd.first <- svd(X)
    V.first <- svd.first$v
  }

  U.fisrt <- X %*% V.first
  X.new <- U.fisrt[,1:out.num]

  if (is.null(V.second)){
    svd.second <- svd(X_sec_mat)
    V.second <- svd.second$v
  }

  U.second <- X_sec_mat %*% V.second
  X.new <- cbind(X.new, U.second[,1:basis.num])
  X.basis <- cbind(X, X_sec_mat)

  for (t in 1:length(X.new[1, ])) {
    X.new[,t] <- (X.new[,t] - mean(X.new[,t])) / (sd(X.new[,t]) + 1e-50)
  }

  return(list(x = X.new, basis = X.basis, V.first = V.first,
              V.second = V.second))
}



########################################################
################ SL with random forest ##################
####### Used for the comparison in Section S5  #########
########################################################


RF_predict <- function(Yt, Xt, St, X.test, S.test, n.tree = 100){
  reg <- Xt
  test.reg <- X.test

  if (is.null(St)){

  }else{
    for (j in 1:(length(unique(St)) - 1)) {
      reg <- cbind(reg, ifelse(St == j, 1, 0))
      test.reg <- cbind(test.reg, ifelse(S.test == j, 1, 0))
    }
  }


  rf_fit <- randomForest(reg, as.factor(Yt), ntree = n.tree)
  ##print(length(test.reg[1,]))
  rf_class <- predict(rf_fit, test.reg)
  rf_class <- as.numeric(rf_class) - 1

  rf_fit_reg <- randomForest(reg, Yt, ntree = n.tree)
  rf_pred <- predict(rf_fit_reg, test.reg)

  return(list(class = rf_class, pred = rf_pred))
}


########################################################
######### Evaluate out-of-sampe AUC, BS, OMR ###########
#### Used for the comparison in Section S5  #############
########################################################

logit_pred <- function(X.test, beta.fit){
  pred_mean <- g.logit(cbind(1, X.test) %*% beta.fit)
  pred_mean[which(is.na(pred_mean))] <- ifelse(cbind(1, X.test) %*% beta.fit > 0, 1, 0)
  return(list(class = ifelse(pred_mean > 0.5, 1, 0), pred = pred_mean))
}

Eva_class <- function(Y.test, results, weight = rep(1 / length(Y.test),
                                                    length(Y.test))){

  auc <- auc(Y.test, results$class)
  mse <- sum((Y.test - results$pred)^2 * weight)
  ae <- sum(abs(Y.test - results$class) * weight)

  return(list(auc = auc, mse = mse, ae = ae))
}


