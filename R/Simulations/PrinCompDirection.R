######################################################################
############ Extract direction of the principle component ############
############ Used for the ML settings in Section S5  #################
######################################################################


extract.pc <- function(X, out.num = 5, basis.num = 10,
                       V.first = NULL, V.second = NULL){
  p <- length(X[1,])
  X_sec_mat <- c()
  for (j in 1:p){
    for (k in (j:p)){
      X_sec_mat <- cbind(X_sec_mat, X[,j] * X[,k])
    }
  }

  if (is.null(V.first)){
    svd.first <- svd(X)
    V.first <- svd.first$v
  }

  U.fisrt <- X %*% V.first
  X.new <- U.fisrt[,1:out.num]

  if (is.null(V.second)){
    svd.second <- svd(X_sec_mat)
    V.second <- svd.second$v
  }

  U.second <- X_sec_mat %*% V.second
  X.new <- cbind(X.new, U.second[,1:basis.num])
  X.basis <- cbind(X, X_sec_mat)

  for (t in 1:length(X.new[1, ])) {
    X.new[,t] <- (X.new[,t] - mean(X.new[,t])) / (sd(X.new[,t]) + 1e-50)
  }

  return(list(x = X.new, basis = X.basis, V.first = V.first,
              V.second = V.second))
}


