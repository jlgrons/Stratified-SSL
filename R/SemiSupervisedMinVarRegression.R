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

SemiSupervisedMinVarRegression <- function(beta_SSL, beta_SL, resids_imp,
                                           resids_beta_SL, X_labeled,
                                           X_unlabeled, epsilon = NULL){

  X_all <- cbind(1, rbind(X_labeled, X_unlabeled))
  X_labeled_int <- cbind(1, X_labeled)

  n_all <- nrow(X_all)
  n_labeled <- nrow(X_labeled_int)
  ones <- c(1, 1)
  p <- ncol(X_labeled)
  weight <- matrix(NA, p + 1, 2)

  if(is.null(epsilon)){epsilon = rep(0, p + 1)}

  # Compute the minimum variance estimator.
  # Note: The residuals have been divided by mean of weights.
  pred_prob_deriv <- c(ExpitDerivative(X_all %*% beta_SSL))
  info_matrix <- solve(t(X_all) %*% (X_all * pred_prob_deriv))
  scaled_info_matrix <- info_matrix * n_all

  imp_IF <- scaled_info_matrix %*% t(X_labeled_int * c(resids_imp))
  beta_SL_IF <- scaled_info_matrix %*% t(X_labeled_int * c(resids_beta_SL))

  for (i in 1:(p + 1)){
    mat <- cbind(imp_IF[i,], beta_SL_IF[i,])
    cov <- solve(t(mat)%*%mat  + epsilon[i]*diag(2))/n_labeled
    w2 <- t(ones) %*% cov %*% ones
    w1 <- t(ones) %*% cov
    weight[i,] <- w1/c(w2)
    if (NA %in% weight[i,]){
      weight[i,] <- c(1, 0)
    }
  }

  # Componentwise minimum variance estimator.
  beta_SSL.w <- beta_SSL*weight[,1] + weight[,2]*beta_SL

  # Standard error of minimum variance estimator.
  A <- solve(t(X_all) %*% (X_all*c(dg.logit(X_all%*%beta_SSL.w))))*n_all;
  w.beta_SSL <- diag(weight[,1]); w.beta_SL =  diag(weight[,2]);
  resids.w <- (w.beta_SSL %*% t(X_labeled_int*c(resids_imp)) +
                w.beta_SL %*% t(X_labeled_int*resids_beta_SL))/n_labeled
  se.beta.w <- sqrt(diag(A %*% (resids.w %*% t(resids.w)) %*% A))

  return(list(beta = beta_SSL.w, weight = weight, se.est = se.beta.w))
}

