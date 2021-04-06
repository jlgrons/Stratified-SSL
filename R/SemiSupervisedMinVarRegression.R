# computes the component-wise minimum variance estimator based on the CV resids
min.var.est = function(beta.ssl, beta.sl, resids.gamma, resids.beta.sl,
                       Xt, Xv, eps = NULL){
  Xv.1 = cbind(1, rbind(Xt,Xv)); Xt.1 = cbind(1,Xt);
  n.v = nrow(Xv.1); n.t = nrow(Xt.1); ones = c(1,1);
  p = ncol(Xt); weight = matrix(NA, p +1, 2);

  if(is.null(eps)){eps = rep(0, p+1)};

  # compute minimum var estimator (note the residuals have been divided by mean of weights)

  A = solve(t(Xv.1) %*% (Xv.1*c(dg.logit(Xv.1%*%beta.ssl))))*n.v;
  T_1 = A %*% t(Xt.1*c(resids.gamma));
  T_2 = A %*% t(Xt.1*c(resids.beta.sl));

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
  beta.ssl.w = beta.ssl*weight[,1] + weight[,2]*beta.sl

  # std error of min var
  A = solve(t(Xv.1) %*% (Xv.1*c(dg.logit(Xv.1%*%beta.ssl.w))))*n.v;
  w.beta.ssl = diag(weight[,1]); w.beta.sl =  diag(weight[,2]);
  resids.w = (w.beta.ssl %*% t(Xt.1*c(resids.gamma)) +
                w.beta.sl %*% t(Xt.1*resids.beta.sl))/n.t
  se.beta.w = sqrt(diag(A %*% (resids.w %*% t(resids.w)) %*% A))


  return(list(beta = beta.ssl.w, weight= weight, se.est = se.beta.w))
}

