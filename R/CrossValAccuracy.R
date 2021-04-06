################################################################################
##### Cross-validated estimator of the model accuracy parameter (BS and OMR) ###
################################################################################

# computes the CV estimates for model evaluation parameters
model.eval.cv = function(basis.x, Xv, Xt, Yt, samp.prob,
                         W, K.fold = 5, rep = 10, c = 0.5, lambda0 = NULL){

  ## basis.x: basis exp matrix
  ## Xt: covariates for labeled
  ## Xv: covariates for unlabeled
  ## Yt: labels
  ## samp.prob: sampling probability
  ## W: weight for the minimum variance estimator
  ## K.fold: number of folds for CV
  ## rep: number of reps for CV
  ## c: cutoff for the overall misclassification rate
  ## lambda0: ridge parameter for the imputation model

  pp = ncol(basis.x); n.t = nrow(Xt); ind.lab = 1:n.t;
  data.all = rbind(Xt,Xv);

  mse.cv.ssl = matrix(NA, nrow = K.fold, ncol = rep);
  ae.cv.ssl = matrix(NA, nrow = K.fold, ncol = rep);

  mse.cv.sl = matrix(NA, nrow = K.fold, ncol = rep);
  ae.cv.sl = matrix(NA, nrow = K.fold, ncol = rep);

  mse.cv.naive = matrix(NA, nrow = K.fold, ncol = rep);
  ae.cv.naive = matrix(NA, nrow = K.fold, ncol = rep);

  mse.cv.dr = matrix(NA, nrow = K.fold, ncol = rep);
  ae.cv.dr = matrix(NA, nrow = K.fold, ncol = rep);


  if(is.null(lambda0)){lambda0 = log(pp)/(floor((K.fold-1)*n.t/K.fold))^1.5};

  for(j in 1:rep){
    set.seed(j)
    ind.cv = split(1:n.t, sample(rep(1:K.fold, floor(n.t/K.fold))));

    for(i in 1:K.fold){

      inds.v = as.vector(unlist(ind.cv[i]));
      inds.t = as.vector(setdiff(as.vector(unlist(ind.cv)), inds.v))

      wg.v = samp.prob[inds.v]; wg.t = samp.prob[inds.t]

      Yt.v = Yt[inds.v]; Xt.v = Xt[inds.v,]; basis.v = basis.x[inds.v, ];
      Yt.t = Yt[inds.t]; Xt.t = Xt[inds.t,]; basis.t = basis.x[inds.t, ]

      beta.tmp = glm.fit.SS(rbind(basis.t, basis.x[-ind.lab, ]),
                            Xt.t, Xv, Yt.t, wg.t, lambda0)

      beta.ssl.t =  beta.tmp$beta.ssl
      beta.sl.t = beta.tmp$beta.sl
      gamma.t =  beta.tmp$gamma
      beta.ssl.w.t = W*beta.ssl.t + (1-W)*beta.sl.t;
      beta.dr.t = beta.tmp$beta.dr

      lp.v = g.logit(cbind(1,Xt.v)%*%beta.ssl.w.t);
      lp.v.ind = ifelse(I(lp.v > c), 1, 0);

      imps.v = cbind(1, basis.v)%*%gamma.t;
      refit.pv.1 <- glm(Yt.v~cbind(lp.v), offset = imps.v,
                        weights = 1/wg.v/mean(1/wg.v), family = 'binomial')$coeff
      refit.pv.2 <-glm(Yt.v~cbind(lp.v.ind), offset = imps.v,
                       weights = 1/wg.v/mean(1/wg.v), family = 'binomial')$coeff

      lp.u = g.logit(cbind(1,data.all)%*%beta.ssl.w.t);
      lp.u.ind = I(lp.u > c);

      imps.pe.1 = g.logit(cbind(1,lp.u)%*% refit.pv.1 + cbind(1, basis.x)%*%gamma.t);
      imps.pe.2 = g.logit(cbind(1,lp.u.ind)%*% refit.pv.2 + cbind(1, basis.x)%*%gamma.t);

      mse.cv.ssl[i,j] = mean(imps.pe.1 + (lp.u - 2*imps.pe.1)*lp.u)
      ae.cv.ssl[i,j] = mean(imps.pe.2 + (lp.u.ind - 2*imps.pe.2)*lp.u.ind)

      lp.v.sl = g.logit(cbind(1,Xt.v) %*% beta.sl.t);
      lp.v.sl.ind = I(lp.v.sl > c);

      mse.cv.sl[i,j] = mean((Yt.v  - lp.v.sl)^2 * 1/wg.v)/mean(1/wg.v)
      ae.cv.sl[i,j] = mean(abs(Yt.v - lp.v.sl.ind) * 1/wg.v)/mean(1/wg.v)

      mse.cv.naive[i,j] = mean((Yt.v  - lp.v.sl)^2)
      ae.cv.naive[i,j] = mean(abs(Yt.v - lp.v.sl.ind))

      lp.v.dr = g.logit(cbind(1, Xt.v) %*% beta.dr.t);
      lp.v.dr.ind = I(lp.v.dr > c);

      mse.cv.dr[i,j] = mean((Yt.v  - lp.v.dr)^2 * 1 / wg.v) / mean(1 / wg.v)
      ae.cv.dr[i,j] = mean(abs(Yt.v - lp.v.dr.ind) * 1 / wg.v) / mean(1 / wg.v)

    }
  }

  return(list(mse.ssl = mean(mse.cv.ssl, na.rm = T),
              ae.ssl = mean(ae.cv.ssl, na.rm = T),
              mse.sl = mean(mse.cv.sl, na.rm = T),
              ae.sl = mean(ae.cv.sl, na.rm = T),
              mse.dr = mean(mse.cv.dr, na.rm = T),
              ae.dr = mean(ae.cv.dr, na.rm = T),
              mse.naive = mean(mse.cv.naive, na.rm = T),
              ae.naive = mean(ae.cv.naive, na.rm = T)))

}
