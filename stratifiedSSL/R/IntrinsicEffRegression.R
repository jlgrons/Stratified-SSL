#' Intrinsic efficiency regression.
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


IntrinsicEffEstBeta <- function(basis_labeled, basis_unlabeled,
                                      X_labeled, X_unlabeled, y,
                                      samp_prob, indx_mom,
                                      lambda0 = NULL, theta_prelim){

  n_labeled <- length(y)
  pp <- ncol(basis)
  p <- ncol(X_labeled)

  if(is.null(lambda0)){lambda0 = log(pp)/n_labeled^1.5}

  dat_all <- rbind(X_labeled, X_unlabeled)
  basis_all <- rbind(basis_labeled, basis_unlabeled)

  A <- crossprod(cbind(1, dat_all),
                 as.vector(dg.logit(cbind(1, dat_all) %*% theta_prelim)) * cbind(1, dat_all))
  A <- A / nrow(dat_all)
  A_inv <- solve(A)

  weights = 1/samp_prob / mean(1/samp_prob)

  d <- ncol(X_labeled)
  theta_est_vec <- rep(0, d + 1)
  gamma_init <- my.ridge(basis_labeled, y, weights = weights, lambda = lambda0)

  for (j in 1:(d + 1)){
    e <- rep(0, d + 1)
    e[j] <- 1
    w <- (cbind(1, X_labeled) %*% A_inv %*% e)^2 * weights^2
    w <- w / mean(w)

    # step 1: basis function regression
    gamma <- my.ridge.weight(basis_labeled, y, w, weights, indx_mom,
                             lambda0 = lambda0,
                             initial = gamma_init)

    # step 2: fit the SS GLM
    imps_basis = g.logit(cbind(1, basis.x) %*% as.vector(gamma))
    dat_all = rbind(X_labeled, X_unlabeled);
    beta.ssl = tryCatch(glm(imps_basis ~ dat_all,
                            family = 'binomial')$coeff, error = function(e) rep(NA, p+1))

    theta_est_vec[j] <- beta.ssl[j]
  }

  return(list(theta = theta_est_vec, gamma = gamma))
}
