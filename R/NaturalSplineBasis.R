# Computes the natural spline basis
ns.basis = function(X, S, nk, basis.type = NULL){

  ## X: covariate matrix
  ## S: the stratum ID
  ## nk: number of knots
  ## basis.type: type of basis (choices under different settings)

  if (basis.type == 'spline'){
    X = as.matrix(X); basis = NULL;

    for(i in 1:ncol(X)){
      X_i = X[,i];

      # quantiles to determine appropriate knots
      knots <- quantile(X_i, seq(0, 1, length = nk));

      # changes quantiles if there aren't enough unique values
      j = 0;
      while(length(unique(knots)) != nk){
        j = j + 1;
        knots <- unique(quantile(X_i, seq(0, 1, length = nk  + j)))
      }

      # compute the natural spline basis
      d_k <- (trunc.cub(X_i, knots[nk-1])) /(knots[nk] - knots[nk-1]);
      evals <- sapply(1:(nk-2), function(ii){
        d_i <- (trunc.cub(X_i, knots[ii])) /(knots[nk] - knots[ii]);
        basis.new <- d_i - d_k})

      # bind together and include original variable
      basis = cbind(basis, cbind(X_i, evals))

    }
    basis.S <- c()
    for (k in 1:(length(unique(S)) - 1)) {
      basis.S <- cbind(basis.S, ifelse(S == (k - 1), 1, 0))
    }
    basis <- cbind(basis, basis.S)
  }

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
      knots <- quantile(X_i, seq(0, 1, length = nk));

      # changes quantiles if there aren't enough unique values
      j = 0;
      while(length(unique(knots)) != nk){
        j = j + 1;
        knots <- unique(quantile(X_i, seq(0, 1, length = nk  + j)))
      }

      # compute the natural spline basis
      d_k <- (trunc.cub(X_i, knots[nk-1])) /(knots[nk] - knots[nk-1]);
      evals <- sapply(1:(nk-2), function(ii){
        d_i <- (trunc.cub(X_i, knots[ii])) /(knots[nk] - knots[ii]);
        basis.new <- d_i - d_k})

      # bind together and include original variable
      basis = cbind(basis, cbind(X_i, evals))

    }
    basis <- cbind(basis, S)
    for (k in 1:length(X[1,])) {
      basis <- cbind(basis, X[,k] * S)
    }

  }

  return(basis)
}

