# Updated: 2021-04-19

#' Component-wise minimum variance semi-supervised regression.
#'
#' @param beta_SSL Semi-supervised regression coefficient vector.
#' @param beta_SL Supervised regression coefficient vector.
#' @param resids_imp Residuals from the imputation model.
#' @param resids_beta_SL Residuals from the supervised regression model.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param epsilon Small offset to help with numerical stability.
#' @export
#' @return List with parameter and SE est as well as estimated weights.
#'

min.var.est = function(beta_SSL, beta_SL, resids_imp, resids_beta_SL,
                       X_labeled, X_unlabeled, epsilon = NULL){

  X_unlabeled.1 = cbind(1, rbind(X_labeled,X_unlabeled))
  X_labeled.1 = cbind(1,X_labeled)

  n.v = nrow(X_unlabeled.1)
  n.t = nrow(X_labeled.1)
  ones = c(1,1)
  p = ncol(X_labeled)
  weight = matrix(NA, p +1, 2)

  if(is.null(epsilon)){epsilon = rep(0, p+1)};

  # compute minimum var estimator (note the residuals have been divided by mean of weights)

  A = solve(t(X_unlabeled.1) %*% (X_unlabeled.1*c(dg.logit(X_unlabeled.1%*%beta_SSL))))*n.v
  T_1 = A %*% t(X_labeled.1*c(resids_imp))
  T_2 = A %*% t(X_labeled.1*c(resids_beta_SL))

  for (i in 1:(p + 1)){
    mat = cbind(T_1[i,], T_2[i,])
    cov = solve(t(mat)%*%mat  + epsilon[i]*diag(2))/n.t
    w2= t(ones) %*% cov %*% ones
    w1 = t(ones) %*% cov
    weight[i,] = w1/c(w2)
    if (NA %in% weight[i,]){
      weight[i,] <- c(1, 0)
    }

  }

  # Componentwise minimum variance estimator.
  beta_SSL.w = beta_SSL*weight[,1] + weight[,2]*beta_SL

  # Standard error of minimum variance estimator.
  A = solve(t(X_unlabeled.1) %*% (X_unlabeled.1*c(dg.logit(X_unlabeled.1%*%beta_SSL.w))))*n.v;
  w.beta_SSL = diag(weight[,1]); w.beta_SL =  diag(weight[,2]);
  resids.w = (w.beta_SSL %*% t(X_labeled.1*c(resids_imp)) +
                w.beta_SL %*% t(X_labeled.1*resids_beta_SL))/n.t
  se.beta.w = sqrt(diag(A %*% (resids.w %*% t(resids.w)) %*% A))


  return(list(beta = beta_SSL.w, weight= weight, se.est = se.beta.w))
}

