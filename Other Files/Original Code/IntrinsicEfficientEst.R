
# For estimating beta:

IntrinsicEfficientEstBeta <- function(basis.x, Xt, Xv, Yt, samp.prob, 
                                      indx_mom, lambda0 = NULL, theta_prelim){
  n.t = length(Yt); pp = ncol(basis.x); p = ncol(Xt); ind.lab = 1:n.t;
  if(is.null(lambda0)){lambda0 = log(pp)/n.t^1.5}
  
  dat.all = rbind(Xt, Xv)
  A <- crossprod(cbind(1, dat.all), 
                 as.vector(dg.logit(cbind(1, dat.all) %*% theta_prelim)) * cbind(1, dat.all))
  A <- A / nrow(dat.all)
  A_inv <- solve(A)
  weights = 1/samp.prob / mean(1/samp.prob)
  
  d <- ncol(Xt)
  theta_est_vec <- rep(0, d + 1)
  gamma_init <- my.ridge(basis.x[ind.lab, ], Yt, weights = weights, lambda = lambda0)
  
  for (j in 1:(d + 1)){
    e <- rep(0, d + 1)
    e[j] <- 1
    w <- (cbind(1, Xt) %*% A_inv %*% e)^2 * weights^2
    w <- w / mean(w)
    
    # step 1: basis function regression
    gamma <- my.ridge.weight(basis.x[ind.lab, ], Yt, w, weights, indx_mom, lambda0 = lambda0,
                             initial = gamma_init)
    
    # step 2: fit the SS GLM
    #res.lab <- Yt - g.logit(cbind(1, basis.x[ind.lab, ]) %*% as.vector(gamma))
    imps.basis = g.logit(cbind(1, basis.x) %*% as.vector(gamma))
    dat.all = rbind(Xt, Xv);
    beta.ssl = tryCatch(glm(imps.basis ~ dat.all, 
                            family = 'binomial')$coeff, error = function(e) rep(NA, p+1))
    
    theta_est_vec[j] <- beta.ssl[j]
  }
  return(list(theta = theta_est_vec, gamma = gamma))
}

  




# For estimating the Brier score:

IntrinsicEfficientEstBS <- function(basis.x, Xt, Xv, Yt, samp.prob, indx_mom,
                                    lambda0 = NULL, theta_prelim){
  dat.all = rbind(Xt, Xv)
  N <- nrow(dat.all)
  n <- nrow(Xt); pp = ncol(basis.x); p = ncol(Xt); ind.lab = 1:n;
  if(is.null(lambda0)){lambda0 = log(pp)/n^1.5}
  
  z_prelim <- as.vector(cbind(1, dat.all) %*% theta_prelim)
  A <- crossprod(cbind(1, dat.all), dg.logit(z_prelim) * cbind(1, dat.all))
  A <- A / N
  A_inv <- solve(A)
  weights = 1 / samp.prob / mean(1 / samp.prob)
  Ddot <- colSums(- 2 * weights * dg.logit(z_prelim[ind.lab]) *
                    (Yt - g.logit(z_prelim[ind.lab])) * cbind(1, Xt)) / n
  
  lp.all_prelim <- g.logit(z_prelim)
  basis.x <- cbind(basis.x, lp.all_prelim)
  indx_mom <- c(indx_mom, ncol(basis.x))
  
  w <- weights^2 * (1 - 2 * lp.all_prelim[ind.lab] + cbind(1, Xt) %*% A_inv %*% Ddot)^2
  w <- as.vector(w / mean(w))
  
  # step 1: basis function regression (for intrinsic estimator)
  gamma_init <- my.ridge(basis.x[ind.lab, ], Yt, weights = weights, lambda = lambda0)
  gamma <- my.ridge.weight(basis.x[ind.lab, ], Yt, w, weights, indx_mom, lambda0 = lambda0,
                           initial = gamma_init)
  gamma <- as.vector(gamma)
  
  # step 2: fit the SS GLM
  imps.basis = g.logit(cbind(1, basis.x) %*% gamma)
  beta.ssl.D = tryCatch(glm(imps.basis ~ dat.all, 
                            family = 'binomial')$coeff, error = function(e) rep(NA, p+1))
  
  # step 3: obtain estimation for D
  lp.all <- g.logit(as.vector(cbind(1, dat.all) %*% beta.ssl.D))
  imps.pe = g.logit(cbind(1, basis.x) %*% gamma);
  mse.ssl = mean(imps.pe + (lp.all - 2*imps.pe)*lp.all)
  return(list(value = mse.ssl, basis = basis.x))
}



# For estimating the OMR:
# h is the bandwidth of the kernel used for estimating the nuisance weights

IntrinsicEfficientEstOMR <- function(basis.x, Xt, Xv, Yt, samp.prob, indx_mom,
                                     lambda0 = NULL, theta_prelim, c = 0.5, h = NULL){
  
  dat.all = rbind(Xt, Xv)
  N <- nrow(dat.all)
  n <- nrow(Xt); pp = ncol(basis.x); p = ncol(Xt); ind.lab = 1:n;
  if(is.null(lambda0)){lambda0 = log(pp)/n^1.5}
  if(is.null(h)){h = 1/n^0.25}
  
  # Use kernel smoothing to estimate the nuisance weights:
  
  z_prelim <- as.vector(cbind(1, dat.all) %*% theta_prelim)
  A <- crossprod(cbind(1, dat.all), dg.logit(z_prelim) * cbind(1, dat.all))
  A <- A / N
  A_inv <- solve(A)
  weights = 1 / samp.prob / mean(1 / samp.prob)
  Kh <- 1 / sqrt(2 * pi) / h * exp(- (g.logit(z_prelim[ind.lab]) - c)^2 / h^2 / 2)
  Ddot <- colSums(weights * dg.logit(z_prelim[ind.lab]) * Kh * (1 - 2 * Yt) * cbind(1, Xt)) / n
  
  lp.all_prelim <- I(g.logit(z_prelim) > c)
  basis.x <- cbind(basis.x, lp.all_prelim)
  indx_mom <- c(indx_mom, ncol(basis.x))
  
  w <- weights^2 * (1 - 2 * lp.all_prelim[ind.lab] + cbind(1, Xt) %*% A_inv %*% Ddot)^2
  w <- as.vector(w / mean(w))
  
  # step 1: basis function regression
  gamma_init <- my.ridge(basis.x[ind.lab, ], Yt, weights = weights, lambda = lambda0)
  gamma <- my.ridge.weight(basis.x[ind.lab, ], Yt, w, weights, indx_mom, lambda0 = lambda0,
                           initial = gamma_init)
  gamma <- as.vector(gamma)
  
  # step 2: fit the SS GLM
  imps.basis = g.logit(cbind(1, basis.x) %*% gamma)
  beta.ssl.D = tryCatch(glm(imps.basis ~ dat.all, 
                            family = 'binomial')$coeff, error = function(e) rep(NA, p+1))
  
  # step 3: obtain estimation for D
  lp.all <- I(g.logit(as.vector(cbind(1, dat.all) %*% beta.ssl.D)) > c)
  imps.pe = g.logit(cbind(1, basis.x) %*% gamma);
  ae.ssl = mean(imps.pe + (lp.all - 2*imps.pe)*lp.all)
  return(list(value = ae.ssl, basis = basis.x))
}





##### Help functions:

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


# Fit penalized and weighted ridge model for intrinsic efficient estimator.

my.ridge.weight <- function(x, y, weights, weights_mom, indx_mom, lambda0 = 1e-04,
                            initial = rep(0, 1 + ncol(x))){
  
  ## x: predictors
  ## y: outcome 
  ## weights: sampling or resampling weights 
  
  if(is.null(weights)){weights = rep(1, length(y))}
  
  # multiple lambdas to make sure convergence happens
  gamma = Newton_glmnet(x, y, weights, weights_mom, indx_mom, lambda0 =lambda0,
                        initial = initial)
  return(gamma)
  
}


# Solve a (regularized) least square problem with expit link under moment constraints.

Newton_glmnet <- function(x, y, weights, weights_mom, indx_mom, 
                          lambda0 = 1e-04, max.iter = 100, tol = 1e-4,
                          initial = rep(0, 1 + ncol(x))){
  error <- Inf
  iter <- 0
  gamma <- initial
  x <- cbind(1, x)
  n <- nrow(x)
  
  sqloss <- mean(weights * (y - g.logit(as.vector(x %*% gamma)))^2) 
  indx_mom <- c(1, 1 + indx_mom)

  while(iter < max.iter & error > tol){
    
    iter <- iter + 1
    gamma_old <- gamma
    sqloss_old <- sqloss
    
    # Update the minimization:
    
    z <- as.vector(x %*% gamma)
    y_ <- y - g.logit(z) + dg.logit(z) * z
    x_ <- as.vector(dg.logit(z)) * x
    xTx <- crossprod(x_, as.vector(weights) * x_) / n + lambda0 / 2 * diag(rep(1, ncol(x_)))
    xTy <- t(x_) %*% (as.vector(y_) * as.vector(weights)) / n
    C <- crossprod(x[,indx_mom], as.vector(weights_mom) * x_) / n
    b <- as.vector(t(x[,indx_mom]) %*% (as.vector(y_) * as.vector(weights_mom)) / n)
    mat_bind <- cbind(rbind(xTx, C), rbind(t(C), matrix(0, nrow(C), nrow(C))))
    vec_bind <- c(xTy, b)
    solution_bind <- solve(mat_bind) %*% vec_bind
    gamma <- solution_bind[1:ncol(x)] 
    sqloss <- mean(weights * (y - g.logit(as.vector(x %*% gamma)))^2)
    
    if (sqloss_old < sqloss){
      gamma <- gamma_old
    }
    error <- sqrt(mean((gamma - gamma_old)^2)) 
  }
  return(gamma)
}


