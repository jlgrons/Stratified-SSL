# Updated: 2021-04-22

#' Standard error estimation based on IF expansion for density ratio method.
#'
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param beta Regression parameter estimate to compute the score function.
#' @param dr_projection Projection from density ratio (DR) method.
#' @param residual Estimated residual.
#' @export
#' @return Standard error estimate.
#'

StdErrorEstimation <- function(X_labeled, X_unlabeled, y, beta, dr_projection,
                               residual){

  X_labeled_intercept <- cbind(1, X_labaled)
  X_all_intercept <- cbind(1, rbind(X_labeled, X_unlabeled))

  n_all <- nrow(X_all_intercept)
  n_labeled <- nrow(X_labeled_intercept)

  dr_residual <- diag(as.vector(residual)) %*% X_labeled_intercept
  score <-  t(dr_residual -
                dr_projection) %*% (dr_residual - dr_projection) / n_labled

  expit_deriv <- c(ExpitDerivative(X_all_intercept %*% beta))
  # Ask Molei about this.
  expit_deriv[which(is.na(expit_deriv))] <- 0

  information <- t(X_all_intercept) %*% (X_all_intercept * expit_deriv) / n_all
  inverse_information <- solve(information)
  variance_est <-  inverse_information %*% score %*% inverse_information

  return(list(std_error = sqrt(diag(inverse_information) / n_labeled),
              inverse_information = inverse_information))
}









se.est.dr <- function(beta, Yt, Xt, Xv, resid, dr_projection){
  Xt.1 = cbind(1, Xt); Xv.1 = cbind(1, rbind(Xt, Xv));
  n.v = nrow(Xv.1); n.t = nrow(Xt.1);

  u.dr <- diag(as.vector(resid)) %*% Xt.1
  S =  t(u.dr - dr_projection) %*% (u.dr - dr_projection) /n.t;

  ddd = c(dg.logit(Xv.1 %*% beta))
  ddd[which(is.na(ddd))] = 0

  I = t(Xv.1) %*% (Xv.1 * ddd)/n.v;

  A = solve(I);
  vars =  A %*% S %*% A;
  return(sqrt(diag(vars)/n.t))

}
