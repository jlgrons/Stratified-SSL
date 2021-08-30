# Updated: 2021-08-30

#' Ridge regression.
#'
#' @param X Covariate matrix.
#' @param y Numeric outcome vector.
#' @param weights Numeric vector of weights.
#' @param weights_mom Numeric vector of mom weights.
#' @param idx_mom Numeric vector of indices.
#' @param lambda0 Penalization parameter.
#' @param initial Numeric vector with an intial value.
#' @export
#' @return Vector containing regression coefficients.
#'

WeightedRidgeRegression <- function(X, y, weights, weights_mom, indx_mom,
                            lambda0 = 1e-04, initial = rep(0, 1 + ncol(x))){

  if(is.null(weights)){weights = rep(1, length(y))}

  # multiple lambdas to make sure convergence happens
  gamma = Newton_glmnet(X, y, weights, weights_mom, indx_mom, lambda0 =lambda0,
                        initial = initial)

  return(gamma)

}
