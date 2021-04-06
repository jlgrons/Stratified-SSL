#######################################################
### SE Est for model accuracy parameters (BS & OMR) ###
#######################################################

# computes the standard error estimates for model evaluation parameters
model.eval.se = function(b, Yt, Xt, Xv, basis.x, samp.prob, W, beta.sl, beta.ssl.w,
                         beta.dr, resids.beta.sl, resids.gamma, resids.beta.dr,
                         A, c = 0.5){

  ## b: number of resamples
  ## Xt: covariates for labeled
  ## Xv: covariates for unlabeled
  ## Yt: labels
  ## basis.x: basis exp matrix
  ## samp.prob: sampling probability
  ## W: weight for the minimum variance estimator
  ## beta.sl: supervised beta
  ## beta.ssl.w: semi-supervised beta (minVar)
  ## resids.beta.sl: CV resids for supervised beta
  ## resids.gamma: CV resids for imputation parameter
  ## A: information matrix based on semi-supervisd beta (minVar)
  ## c: cutoff for the overall misclassification rate

  wgt = sapply(1:b, function(kk) 4*rbeta(n.t, 0.5, 1.5))

  resids.beta.sl.p = ((wgt-1)*resids.beta.sl)
  resids.gamma.p = ((wgt-1)*resids.gamma)
  resids.beta.dr.p = (wgt-1)*resids.beta.dr
  proj.dr.p = t(proj.dr) %*% (wgt-1)

  Xt.1 = cbind(1, Xt)
  T_1.p = A %*% t(Xt.1)%*%resids.gamma.p/n.t;
  T_2.p = A %*% t(Xt.1)%*%resids.beta.sl.p/n.t;

  T_1.dr.p = A %*% t(Xt.1) %*% resids.beta.dr.p / n.t
  T_2.dr.p = A %*% proj.dr.p / n.t

  beta.ssl.w.p = beta.ssl.w + diag(W[,1]) %*% T_1.p + diag(W[,2]) %*% T_2.p;
  beta.sl.p = beta.sl + T_2.p
  beta.dr.p = beta.dr + T_1.dr.p - T_2.dr.p

  ind.lab = 1:n.t;

  gamma.p = sapply(1:b, function(kk){
    my.ridge(basis.x[ind.lab, ], Yt, weights = wgt[,kk]/samp.prob/mean(wgt[,kk]/samp.prob),
             lambda = log(ncol(basis.x))/(n.t^1.5))
  })

  perts = lapply(1:b, function(kk){model.eval.ap(Yt, beta.sl.p[,kk], beta.ssl.w.p[,kk], gamma.p[,kk],
                                                 beta.dr.p[,kk], Xt, Xv, basis.x, samp.prob, V = wgt[,kk])})

  ssl.pert.mse = sapply(perts, "[[", 1);
  ssl.pert.ae = sapply(perts, "[[", 2);

  sl.pert.mse = sapply(perts, "[[", 3);
  sl.pert.ae = sapply(perts, "[[", 4);

  dr.pert.mse = sapply(perts, "[[", 5);
  dr.pert.ae = sapply(perts, "[[", 6);

  return(list(ssl.pert.mse = ssl.pert.mse, sl.pert.mse = sl.pert.mse,
              dr.pert.mse = dr.pert.mse, ssl.pert.ae = ssl.pert.ae,
              sl.pert.ae = sl.pert.ae, dr.pert.ae = dr.pert.ae,
              gamma.p = gamma.p, beta.ssl.w.p = beta.ssl.w.p,
              beta.sl.p = beta.sl.p, beta.dr.p = beta.dr.p))

}
