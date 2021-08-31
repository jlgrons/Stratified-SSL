
DensityRatioAccuracy <- function(Yt, beta.dr, Xt, 
                                 samp.prob, DR.est = NULL, V = NULL, c = 0.5){
  
  if(is.null(DR.est)){DR.est = rep(1, length(Yt))}

  weight = 1 / samp.prob / mean(1 / samp.prob);
  lp.t.dr = Expit(cbind(1, Xt) %*% beta.dr);
  lp.t.dr.ind = I(lp.t.dr > c);
  
  mse.dr = mean((Yt - lp.t.dr)^2 * DR.est * weight)
  ae.dr = mean(abs(Yt - lp.t.dr.ind) * DR.est * weight)
  
  return(list(mse.dr = mse.dr, ae.dr = ae.dr))
}