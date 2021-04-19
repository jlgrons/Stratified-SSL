# Updated: 2021-04-19

#' Component-wise minimum variance semi-supervised regression.
#'
#' @param beta_SSL Semi-supervised
#' @param beta_SL Basis matrix for unlabeled data set.
#' @param resids_imp Basis matrix for labeled data set.
#' @param resids_beta_SL Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of weights.
#' @param lambda Penalization parameter for initial ridge estimator.
#' @export
#' @return Vector containing regression coefficients.
#'

min.var.est = function(beta_SSL, beta_SL, resids_imp, resids_beta_SL,
                       X_labeled, X_unlabeled, eps = NULL){

  X_unlabeled.1 = cbind(1, rbind(X_labeled,X_unlabeled)); X_labeled.1 = cbind(1,X_labeled);
  n.v = nrow(X_unlabeled.1); n.t = nrow(X_labeled.1); ones = c(1,1);
  p = ncol(X_labeled); weight = matrix(NA, p +1, 2);

  if(is.null(eps)){eps = rep(0, p+1)};

  # compute minimum var estimator (note the residuals have been divided by mean of weights)

  A = solve(t(X_unlabeled.1) %*% (X_unlabeled.1*c(dg.logit(X_unlabeled.1%*%beta_SSL))))*n.v;
  T_1 = A %*% t(X_labeled.1*c(resids_imp));
  T_2 = A %*% t(X_labeled.1*c(resids_beta_SL));

  for (i in 1:(p + 1)){
    mat = cbind(T_1[i,], T_2[i,]);
    cov = solve(t(mat)%*%mat  + eps[i]*diag(2))/n.t;
    w2= t(ones) %*% cov %*% ones
    w1 = t(ones) %*% cov
    weight[i,] = w1/c(w2)
    if (NA %in% weight[i,]){
      weight[i,] <- c(1, 0)
    }

  }

  # min var estimator
  beta_SSL.w = beta_SSL*weight[,1] + weight[,2]*beta_SL

  # std error of min var
  A = solve(t(X_unlabeled.1) %*% (X_unlabeled.1*c(dg.logit(X_unlabeled.1%*%beta_SSL.w))))*n.v;
  w.beta_SSL = diag(weight[,1]); w.beta_SL =  diag(weight[,2]);
  resids.w = (w.beta_SSL %*% t(X_labeled.1*c(resids_imp)) +
                w.beta_SL %*% t(X_labeled.1*resids_beta_SL))/n.t
  se.beta.w = sqrt(diag(A %*% (resids.w %*% t(resids.w)) %*% A))


  return(list(beta = beta_SSL.w, weight= weight, se.est = se.beta.w))
}

