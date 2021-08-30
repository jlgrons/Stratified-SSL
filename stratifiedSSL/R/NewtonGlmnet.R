# Updated: 2021-08-30

#' Solve a (regularized) least square problem
#' with expit link under moment constraints.
#'
#' @param X Covariate matrix.
#' @param y Numeric outcome vector.
#' @param weights Numeric vector of weights.
#' @param weights_mom Numeric vector of mom weights.
#' @param idx_mom Numeric vector of indices.
#' @param lambda0 Penalization parameter.
#' @param max_iter Max number of iterations.
#' @param tol Desired tolerance level.
#' @param initial Numeric vector with an intial value.
#' @export
#' @return Vector containing regression coefficients.
#'

NewtonGlmnet <- function(x, y, weights, weights_mom, indx_mom,
                          lambda0 = 1e-04, max_iter = 100, tol = 1e-4,
                          initial = rep(0, 1 + ncol(x))){
  error <- Inf
  iter <- 0
  gamma <- initial
  x <- cbind(1, X)
  n <- nrow(X)

  sqloss <- mean(weights * (y - g.logit(as.vector(X %*% gamma)))^2)
  indx_mom <- c(1, 1 + indx_mom)

  while(iter < max_iter & error > tol){

    iter <- iter + 1
    gamma_old <- gamma
    sqloss_old <- sqloss

    # Update the minimization:

    z <- as.vector(X %*% gamma)
    y_ <- y - g.logit(z) + dg.logit(z) * z
    x_ <- as.vector(dg.logit(z)) * X
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
