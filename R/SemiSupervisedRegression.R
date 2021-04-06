# fits the glm using unlabeled and labeled data & computes corresponding supervised estimator
glm.fit.SS = function(basis.x, Xt, Xv, Yt, samp.prob, lambda0 = NULL){

  ## basis.x: basis matrix with first n.t rows corr to labeled data and unlabeled rows n.t+1 to n.t+n.v
  ## Xt: covariates for labeled
  ## Xt: covariates for labeled
  ## Yt: labels
  ## samp.prob: vector of sampling weights
  ## lambda0: ridge parameter for SS estimator

  n.t = length(Yt); pp = ncol(basis.x); p = ncol(Xt); ind.lab = 1:n.t;

  if(is.null(lambda0)){lambda0 = log(pp)/n.t^1.5};

  weights = 1/samp.prob/mean(1/samp.prob)

  # step 1: basis function regression
  gamma = my.ridge(basis.x[ind.lab, ], Yt, weights = weights, lambda = lambda0);

  # step 2: fit the SS GLM
  imps.basis = g.logit(cbind(1, basis.x) %*% gamma);
  dat.all = rbind(Xt, Xv);
  beta.ssl = tryCatch(glm(imps.basis ~ dat.all,
                          family = 'binomial')$coeff, error = function(e) rep(NA, p+1))

  # fit supervised estimate for comparison
  beta.sl = tryCatch(glm(Yt ~ Xt, family = 'binomial',
                         weights = weights)$coeff, error = function(e) rep(NA, p+1));

  beta.naive = tryCatch(glm(Yt ~ Xt, family = 'binomial')$coeff,
                        error = function(e) rep(NA, p+1));

  # Fit the density ratio estimator

  n.v <- length(basis.x[,1])
  basis.phi <- cbind(rep(1, n.v), basis.x)
  phi.t <- basis.phi[1:n.t, ]
  phi.v <- basis.phi[- c(1:n.t), ]
  dim.basis <- length(basis.phi[1, ])

  phiT.phi <- t(basis.phi) %*% basis.phi / (n.t + n.v)
  E.phi <- t(phi.t) %*% weights / n.t
  theta.ratio <- solve(phiT.phi + diag(rep(0.01, ncol(phiT.phi)))) %*% (rowMeans(t(phi.v)) - E.phi)
  density.ratio <- exp(phi.t %*% theta.ratio)
  weights.dr <- density.ratio * weights

  beta.dr <- tryCatch(glm(Yt~Xt, family = 'binomial',
                          weights = weights.dr)$coeff, error = function(e) rep(NA, p+1));

  u.dr <- diag(as.vector(Yt - g.logit(cbind(1, Xt) %*% beta.dr))) %*% cbind(1, Xt)
  E.uT.phi <- t(u.dr) %*% diag(as.vector(weights)) %*% phi.t / sum(weights)
  proj.coef.dr <- E.uT.phi %*% solve(phiT.phi + diag(rep(0.01, ncol(phiT.phi))))
  proj.dr <- diag(as.vector(weights)) %*% phi.t %*% t(proj.coef.dr)

  return(list(beta.ssl = beta.ssl, beta.sl = beta.sl, gamma = gamma,
              beta.dr = beta.dr, proj.dr = proj.dr,
              beta.naive = beta.naive))
}
