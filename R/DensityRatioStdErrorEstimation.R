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
