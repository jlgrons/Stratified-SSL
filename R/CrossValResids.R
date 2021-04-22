# Updated: 2021-04-22

#' Apparent estimates for brier score (MSE) and misclassification rate (OMR).
#'
#' @param basis_labeled Basis matrix for labeled data set.
#' @param bais_unlabeled Basis matrix for unlabeled data set.
#' @param X_labeled Covariate matrix for labeled data set.
#' @param X_unlabeled Covariate matrix for unlabeled data set.
#' @param y Numeric outcome vector.
#' @param samp_prob Numeric vector of weights.
#' @param num_folds Number of folds for cross-validation.
#' @param lambda Regularization parameter for the imputation model.
#' @export
#' @return Cross-validated residuals.
#'

CrossValResids <- function(basis_labeled, basis_unlabeled, X_labeled,
                           X_unlabeled, y, samp_prob, num_folds, lambda = NULL){

  n_labeled <- nrow(X_labeled)
  p_basis <- ncol(basis_labeled)
  ind_cv <- split(1:n_labeled, sample(rep(1:num_folds,
                                          floor(n_labeled / num_folds))))

  if(is.null(lambda)){
    lambda <- log(p_basis) / (floor((num_folds - 1) * n_labeled/num_folds))^1.5
    }

  # K estimates of supervised and SS beta & gamma with kth fold removed
  cv.ests <- lapply(1:num_folds, function(kk) {
    inds_fold <- as.vector(unlist(ind_cv[-kk]))
    y_fold <- y[inds_fold]
    samp_prob_fold <- samp_prob[inds_fold]
    SemiSupervisedRegression(basis_labeled[inds_fold,], basis_unlabeled,
                             X_labeled[inds_fold,], X_unlabeled, y_fold,
                             samp_prob_fold, lambda = lambda)
  })


  beta.ssl.cv = sapply(cv.ests, "[[", 1);
  beta.sl.cv = sapply(cv.ests, "[[", 2);
  gamma.cv = sapply(cv.ests, "[[", 3);
  beta.dr.cv = sapply(cv.ests, "[[", 4);

  # Basis corresponding to the labeled set
  basis.lab = basis.x[ind.lab, ];

  # K sets of residuals based on the kth fold and the beta with the kth fold removed
  res.cv = lapply(1:num_folds,function(kk){inds.v = as.vector(unlist(ind_cv[kk]));
  beta.ssl.cv.tmp = beta.ssl.cv[,kk]; beta.sl.cv.tmp = beta.sl.cv[,kk];
  gamma.cv.tmp = gamma.cv[,kk]; beta.dr.cv.tmp = beta.dr.cv[,kk];
  samp.prob.v = samp.prob[inds.v];
  Yt.v = Yt[inds.v]; X_labeled.v = X_labeled[inds.v,]; basis.v = basis.lab[inds.v, ];
  pred.b.ssl = g.logit(cbind(1, X_labeled.v) %*% beta.ssl.cv.tmp);
  pred.b.sl = g.logit(cbind(1, X_labeled.v) %*% beta.sl.cv.tmp);
  pred.gamma = g.logit(cbind(1, basis.v) %*% gamma.cv.tmp);
  pred.b.dr = g.logit(cbind(1, X_labeled.v) %*% beta.dr.cv.tmp);
  resids = list(beta.ssl = (Yt.v - pred.b.ssl)*(1/samp.prob.v),
                beta.sl = (Yt.v - pred.b.sl)*(1/samp.prob.v),
                gamma = (Yt.v - pred.gamma)*(1/samp.prob.v),
                beta.dr = (Yt.v - pred.b.dr)*(1/samp.prob.v))
  })

  inds.ord = order(unlist(ind_cv)); wgt = mean(1/samp.prob);

  # Reorder the residuals and divide by the mean of the weights
  resids.beta.ssl = c(unlist(sapply(res.cv, "[[", 1)))[inds.ord]/wgt;
  resids.beta.sl = c(unlist(sapply(res.cv, "[[", 2)))[inds.ord]/wgt;
  resids.gamma = c(unlist(sapply(res.cv, "[[", 3)))[inds.ord]/wgt;
  resids.beta.dr = c(unlist(sapply(res.cv, "[[", 4)))[inds.ord]/wgt;


  return(list(resids.beta.ssl = resids.beta.ssl, resids.beta.sl = resids.beta.sl,
              resids.gamma = resids.gamma, resids.beta.dr = resids.beta.dr,
              beta.ssl.cv = beta.ssl.cv, beta.sl.cv = beta.sl.cv,
              gamma.cv = gamma.cv, beta.dr.cv = beta.dr.cv))
}



