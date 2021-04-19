# Updated: 2021-04-19

#' Density ratio regression_
#'
#' @param basis_labeled Basis matrix for labeled data set_
#' @param bais_unlabeled Basis matrix for unlabeled data set_
#' @param X_labeled Covariate matrix for labeled data set_
#' @param X_unlabeled Covariate matrix for unlabeled data set_
#' @param y Numeric outcome vector_
#' @param samp_prob Numeric vector of weights_
#' @param lambda Penalization parameter for initial ridge estimator_
#' @export
#' @return Vector containing regression coefficients_
#'

DensityRatioRegression <- function(basis_labeled, basis_unlabeled, X_labeled,
                                     X_unlabeled, y, samp_prob, lambda = NULL){

  num_labeled <- nrow(basis_labeled)
  basis_all <- rbind(basis_labeled, basis_unlabeled)
  n_all <- length(basis_x[,1])

  basis_phi <- cbind(rep(1, n_all), basis_all)
  phi_labeled <- basis_phi[1:num_labeled, ]
  phi_unlabeled <- basis_phi[-c(1:num_labeled), ]
  dim_basis <- length(basis_phi[1, ])

  phiT_phi <- t(basis_phi) %*% basis_phi / (num_labeled + n_all)
  E_phi <- t(phi_labeled) %*% weights / num_labeled
  theta_ratio <- solve(phiT_phi + diag(rep(0.01, ncol(phiT_phi)))) %*% (rowMeans(t(phi_unlabeled)) - E_phi)
  density_ratio <- exp(phi_labeled %*% theta_ratio)
  weights_dr <- density_ratio * weights

  beta_dr <- tryCatch(glm(y~X_labeled, family = 'binomial',
                          weights = weights_dr)$coeff, error = function(e) rep(NA, p+1));

  u_dr <- diag(as_vector(y - g_logit(cbind(1, X_labeled) %*% beta_dr))) %*% cbind(1, X_labeled)
  E_uT_phi <- t(u_dr) %*% diag(as_vector(weights)) %*% phi_labeled / sum(weights)
  proj_coef_dr <- E_uT_phi %*% solve(phiT_phi + diag(rep(0.01, ncol(phiT_phi))))
  proj_dr <- diag(as_vector(weights)) %*% phi_labeled %*% t(proj_coef_dr)


  return(list(beta_dr = beta_dr, proj_dr = proj_dr))
}
