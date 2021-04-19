# Updated: 2021-04-19

#' Computes the natural spline basis.
#'
#' @param X Covariate matrix.
#' @param S Vector with stratum id.
#' @param num_knots Number of knots.
#'@param num_knots Number of knots.
#' @export
#' @return Matrix containing basis.
#'

AltenativeBasis <- function(X, S, num_knots){

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
      knots <- quantile(X_i, seq(0, 1, length = num_knots));

      # changes quantiles if there aren't enough unique values
      j = 0;
      while(length(unique(knots)) != num_knots){
        j = j + 1;
        knots <- unique(quantile(X_i, seq(0, 1, length = num_knots  + j)))
      }

      # compute the natural spline basis
      d_k <- (trunc.cub(X_i, knots[num_knots-1])) /(knots[num_knots] - knots[num_knots-1]);
      evals <- sapply(1:(num_knots-2), function(ii){
        d_i <- (trunc.cub(X_i, knots[ii])) /(knots[num_knots] - knots[ii]);
        basis.new <- d_i - d_k})

      # bind together and include original variable
      basis = cbind(basis, cbind(X_i, evals))

    }

    basis <- cbind(basis, S)
    for (k in 1:length(X[1,])) {
      basis <- cbind(basis, X[,k] * S)
    }
  }



}
