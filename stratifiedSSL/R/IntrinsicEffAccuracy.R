#' Intrinsic efficiency OMR.
#'
#' @param basis_labeled Basis matrix for labeled data set.
#' @param bais_unlabeled Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of weights.
#' @param indx_mom Numeric vector of indices.
#' @param lambda0 Penalization parameter for initial ridge estimator.
#' @param theta_prelim Initial estimate.
#' @param threshold Threshold for OMR.
#' @param h Bandwidth of the kernel for nuisance weight estimation.
#' @export
#' @return Vector containing regression coefficients.
#'

IntrinsicEfficientEstOMR <- function(basis_labeled, basis_unlabeled,
                                     X_labeled, X_unlabeled, y, samp_prob,
                                     indx_mom, lambda0 = NULL, theta_prelim,
                                     threshold = 0.5, h = NULL){

  dat_all <- rbind(X_labeled, X_unlabeled)
  basis_all <- rbind(basis_labeled, basis_unlabeled)
  N <- nrow(dat_all)
  n <- nrow(X_labeled)
  pp <- ncol(basis_labeled)
  p <- ncol(X_labeled)
  ind.lab <- 1:n

  if(is.null(lambda0)){lambda0 = log(pp)/n^1.5}
  if(is.null(h)){h = 1/n^0.25}

  # Use kernel smoothing to estimate the nuisance weights:

  z_prelim <- as.vector(cbind(1, dat_all) %*% theta_prelim)
  A <- crossprod(cbind(1, dat_all), dg.logit(z_prelim) * cbind(1, dat_all))
  A <- A / N
  A_inv <- solve(A)
  weights = 1 / samp_prob / mean(1 / samp_prob)
  Kh <- 1 / sqrt(2 * pi) / h * exp(- (Expit(z_prelim[ind.lab]) - c)^2 / h^2 / 2)
  Ddot <- colSums(weights * dg.logit(z_prelim[ind.lab]) * Kh * (1 - 2 * y) * cbind(1, X_labeled)) / n

  lp.all_prelim <- I(Expit(z_prelim) > threshold)
  basis_all <- cbind(basis_all, lp.all_prelim)
  indx_mom <- c(indx_mom, ncol(basis_all))

  w <- weights^2 * (1 - 2 * lp.all_prelim[ind.lab] + cbind(1, X_labeled) %*% A_inv %*% Ddot)^2
  w <- as.vector(w / mean(w))

  # step 1: basis function regression
  gamma_init <- my.ridge(basis_all[ind.lab, ], y, weights = weights, lambda = lambda0)
  gamma <- my.ridge.weight(basis_all[ind.lab, ], y, w, weights, indx_mom, lambda0 = lambda0,
                           initial = gamma_init)
  gamma <- as.vector(gamma)

  # step 2: fit the SS GLM
  imps.basis = Expit(cbind(1, basis_all) %*% gamma)
  beta.ssl.D = tryCatch(glm(imps.basis ~ dat_all,
                            family = 'binomial')$coeff, error = function(e) rep(NA, p+1))

  # step 3: obtain estimation for D
  lp.all <- I(Expit(as.vector(cbind(1, dat_all) %*% beta.ssl.D)) > threshold)
  imps.pe = Expit(cbind(1, basis_all) %*% gamma);
  ae.ssl = mean(imps.pe + (lp.all - 2*imps.pe)*lp.all)
  return(list(value = ae.ssl, basis = basis_all))
}


#' Intrinsic efficiency BS.
#'
#' @param basis_labeled Basis matrix for labeled data set.
#' @param bais_unlabeled Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of weights.
#' @param indx_mom Numeric vector of indices.
#' @param lambda0 Penalization parameter for initial ridge estimator.
#' @param theta_prelim Initial estimate.
#' @export
#' @return Vector containing regression coefficients.
#'

IntrinsicEfficientEstBS <- function(basis_labeled, basis_unlabeled,
                                    X_labeled, X_unlabeled, y,
                                    samp_prob, indx_mom,
                                    lambda0 = NULL, theta_prelim){

  dat_all <- rbind(X_labeled, X_unlabeled)
  basis_all <- rbind(basis_labeled, basis_unlabeled)
  N <- nrow(dat_all)
  n <- nrow(X_labeled)
  pp <- ncol(basis_labeled)
  p <- ncol(X_labeled)

  ind.lab = 1:n;
  if(is.null(lambda0)){lambda0 = log(pp)/n^1.5}

  z_prelim <- as.vector(cbind(1, dat_all) %*% theta_prelim)
  A <- crossprod(cbind(1, dat_all), dg.logit(z_prelim) * cbind(1, dat_all))
  A <- A / N
  A_inv <- solve(A)
  weights = 1 / samp_prob / mean(1 / samp_prob)
  Ddot <- colSums(- 2 * weights * dg.logit(z_prelim[ind.lab]) *
                    (y - Expit(z_prelim[ind.lab])) * cbind(1, X_labeled)) / n

  lp.all_prelim <- Expit(z_prelim)
  basis_all <- cbind(basis_all, lp.all_prelim)
  indx_mom <- c(indx_mom, ncol(basis_all))

  w <- weights^2 * (1 - 2 * lp.all_prelim[ind.lab] + cbind(1, X_labeled) %*% A_inv %*% Ddot)^2
  w <- as.vector(w / mean(w))

  # step 1: basis function regression (for intrinsic estimator)
  gamma_init <- my.ridge(basis_all[ind.lab, ], y, weights = weights, lambda = lambda0)
  gamma <- my.ridge.weight(basis_all[ind.lab, ], y, w, weights, indx_mom, lambda0 = lambda0,
                           initial = gamma_init)
  gamma <- as.vector(gamma)

  # step 2: fit the SS GLM
  imps.basis = Expit(cbind(1, basis_all) %*% gamma)
  beta.ssl.D = tryCatch(glm(imps.basis ~ dat_all,
                            family = 'binomial')$coeff, error = function(e) rep(NA, p+1))

  # step 3: obtain estimation for D
  lp.all <- Expit(as.vector(cbind(1, dat_all) %*% beta.ssl.D))
  imps.pe = Expit(cbind(1, basis_all) %*% gamma);
  mse.ssl = mean(imps.pe + (lp.all - 2*imps.pe)*lp.all)

  return(list(value = mse.ssl, basis = basis_all))
}
