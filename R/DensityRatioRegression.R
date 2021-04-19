# Updated: 2021-04-19

#' Density ratio regression.
#'
#' @param basis_labeled Basis matrix for labeled data set.
#' @param bais_unlabeled Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of weights.
#' @param lambda Penalization parameter for initial ridge estimator.
#' @export
#' @return Vector containing regression coefficients.
#'

DensityRatioRegression <- function(basis_labeled, basis_unlabeled, X_labeled,
                                     X_unlabeled, y, samp_prob, lambda = NULL){
  basis_all <- rbind(basis_labeled, basis_unlabeled)
  n_all <- length(basis.x[,1])
  basis.phi <- cbind(rep(1, n_all), basis.x)
  phi.t <- basis.phi[1:num_labeled, ]
  phi.v <- basis.phi[- c(1:num_labeled), ]
  dim.basis <- length(basis.phi[1, ])

  phiT.phi <- t(basis.phi) %*% basis.phi / (num_labeled + n_all)
  E.phi <- t(phi.t) %*% weights / num_labeled
  theta.ratio <- solve(phiT.phi + diag(rep(0.01, ncol(phiT.phi)))) %*% (rowMeans(t(phi.v)) - E.phi)
  density.ratio <- exp(phi.t %*% theta.ratio)
  weights.dr <- density.ratio * weights

  beta.dr <- tryCatch(glm(y~X_labeled, family = 'binomial',
                          weights = weights.dr)$coeff, error = function(e) rep(NA, p+1));

  u.dr <- diag(as.vector(y - g.logit(cbind(1, X_labeled) %*% beta.dr))) %*% cbind(1, X_labeled)
  E.uT.phi <- t(u.dr) %*% diag(as.vector(weights)) %*% phi.t / sum(weights)
  proj.coef.dr <- E.uT.phi %*% solve(phiT.phi + diag(rep(0.01, ncol(phiT.phi))))
  proj.dr <- diag(as.vector(weights)) %*% phi.t %*% t(proj.coef.dr)


  return(list(beta_SSL = beta_SSL, beta_SL = beta_SL,
              beta_SL_unweighted = beta_SL_unweighted,
              beta_imp = gamma,
              beta.dr = beta.dr, proj.dr = proj.dr,
              beta.naive = beta.naive))
}
