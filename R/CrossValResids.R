#need to change to dg.logit!

# Computes cross-validated (CV) residuals based on supervised and SS ests of the regression parameter
# as well as the regression parameter from the imputation model
resids.cv = function(basis.x, Xt, Xv, Yt, samp.prob, K.fold, lambda0 = NULL){

  ## basis.x: basis matrix with first n.t rows corr to labeled data and unlabeled rows n.t+1 to n.t+n.v
  ## Xt: covariates for labeled
  ## Xv: covariates for unlabeled
  ## Yt: labels
  ## samp.prob: vector of sampling weights
  ## K.fold: # of folds for CV
  ## lambda0: ridge parameter for SS estimator

  n.t = nrow(Xt); ind.cv = split(1:n.t, sample(rep(1:K.fold, floor(n.t/K.fold))));
  qq = length(ind.cv); ind.lab = 1:n.t; pp = ncol(basis.x);

  if(is.null(lambda0)){lambda0 = log(pp)/(floor((K.fold-1)*n.t/K.fold))^1.5};

  # K estimates of supervised and SS beta & gamma with kth fold removed
  cv.ests = lapply(1:qq, function(kk) {inds.t = as.vector(unlist(ind.cv[-kk]));
  Yt.t = Yt[inds.t]; samp.prob.t = samp.prob[inds.t];
  basis.cv = rbind(basis.x[inds.t,], basis.x[-ind.lab,]);
  glm.fit.SS(basis.cv,  Xt[inds.t,], Xv, Yt.t, samp.prob.t, lambda0 = lambda0)});

  beta.ssl.cv = sapply(cv.ests, "[[", 1);
  beta.sl.cv = sapply(cv.ests, "[[", 2);
  gamma.cv = sapply(cv.ests, "[[", 3);
  beta.dr.cv = sapply(cv.ests, "[[", 4);

  # Basis corresponding to the labeled set
  basis.lab = basis.x[ind.lab, ];

  # K sets of residuals based on the kth fold and the beta with the kth fold removed
  res.cv = lapply(1:qq,function(kk){inds.v = as.vector(unlist(ind.cv[kk]));
  beta.ssl.cv.tmp = beta.ssl.cv[,kk]; beta.sl.cv.tmp = beta.sl.cv[,kk];
  gamma.cv.tmp = gamma.cv[,kk]; beta.dr.cv.tmp = beta.dr.cv[,kk];
  samp.prob.v = samp.prob[inds.v];
  Yt.v = Yt[inds.v]; Xt.v = Xt[inds.v,]; basis.v = basis.lab[inds.v, ];
  pred.b.ssl = g.logit(cbind(1, Xt.v) %*% beta.ssl.cv.tmp);
  pred.b.sl = g.logit(cbind(1, Xt.v) %*% beta.sl.cv.tmp);
  pred.gamma = g.logit(cbind(1, basis.v) %*% gamma.cv.tmp);
  pred.b.dr = g.logit(cbind(1, Xt.v) %*% beta.dr.cv.tmp);
  resids = list(beta.ssl = (Yt.v - pred.b.ssl)*(1/samp.prob.v),
                beta.sl = (Yt.v - pred.b.sl)*(1/samp.prob.v),
                gamma = (Yt.v - pred.gamma)*(1/samp.prob.v),
                beta.dr = (Yt.v - pred.b.dr)*(1/samp.prob.v))
  })

  inds.ord = order(unlist(ind.cv)); wgt = mean(1/samp.prob);

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

# computes a standard error estimate using the influence function of
# solution to a simple estimating equation
se.est = function(bhat, Yt, Xt, Xv, resid){

  ## bhat: parameter estimate
  ## Xt: covariates for labeled
  ## Xv: covariates for unlabeled
  ## Yt: labels
  ## resids: vector of CV/non-CV residuals to compute score

  Xt.1 = cbind(1,Xt); Xv.1 = cbind(1, rbind(Xt,Xv));
  n.v = nrow(Xv.1); n.t = nrow(Xt.1);

  S =  t(Xt.1) %*% (Xt.1 * c(resid)^2)/n.t;
  ddd = c(dg.logit(Xv.1%*%bhat))
  ddd[which(is.na(ddd))] = 0
  I = t(Xv.1) %*% (Xv.1 * ddd)/n.v;

  A = solve(I);
  vars =  A %*% S %*% A;
  return(list(se = sqrt(diag(vars)/n.t), A = A))
}


se.est.dr <- function(bhat, Yt, Xt, Xv, resid, proj.dr){
  Xt.1 = cbind(1, Xt); Xv.1 = cbind(1, rbind(Xt, Xv));
  n.v = nrow(Xv.1); n.t = nrow(Xt.1);

  u.dr <- diag(as.vector(resid)) %*% Xt.1
  S =  t(u.dr - proj.dr) %*% (u.dr - proj.dr) /n.t;

  ddd = c(dg.logit(Xv.1 %*% bhat))
  ddd[which(is.na(ddd))] = 0

  I = t(Xv.1) %*% (Xv.1 * ddd)/n.v;

  A = solve(I);
  vars =  A %*% S %*% A;
  return(sqrt(diag(vars)/n.t))

}

