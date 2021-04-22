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



