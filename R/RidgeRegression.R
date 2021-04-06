# fits the penalized ridge model for semi-supervised estimate
my.ridge = function(x, y, weights = NULL, lambda0 = 1e-04){

  ## x: predictors
  ## y: outcome
  ## weights: sampling or resampling weights

  if(is.null(weights)){weights = rep(1, length(y))}

  # multiple lambdas to make sure convergence happens
  gamma = tryCatch((coef(glmnet(x, y, weights = weights, alpha = 0,
                                lambda = seq(lambda0, 100*lambda0, length.out = 100),
                                family = 'binomial'))),
                   error = function(e) rep(NA, ncol(x) + 1 ))

  # take the value at lambda0
  gamma = gamma[, ncol(gamma)]
  return(gamma)
}
