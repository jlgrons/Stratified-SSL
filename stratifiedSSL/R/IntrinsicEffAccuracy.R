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
