set.seed(92047)
setwd('~/Desktop/Stratified-SSL/R/Original Code')
source('function.R')
source('function_response.R')
# Outcome model is wrong and imputation model is correct (Section 7 (ii))

beta.true.all <- c()
beta.sl.all <- c()
beta.ssl.all <- c()
beta.ssl.2.all <- c()

for(i in 1:num_sims){
  data <- data_gen(n.t, N, p, rho, model.spec = 'IC', strata_num = strata.num)

  # Generate the true regression parameters
  data.v <- data_gen(n.t, 100000-n.t, p, rho, model.spec = 'IC')
  X0.v <- data.v$X0
  Y.v <- data.v$Y
  beta.true = glm(Y.v~X0.v, family = 'binomial')$coeff


  # Labled Data + Weights

  Xt <- data$Xt
  Yt <- data$Yt
  Xv <- data$Xv
  samp.prob <- data$samp.prob
  Strata_ID <- data$S
  St <- data$St

  # Basis expansion
  knots = 3;
  dat.all = rbind(Xt, Xv)
  basis.x = ns.basis(dat.all, Strata_ID, knots, basis.type = basis.type)


  # Supervised and SS estimates of regression parameter
  glm.fit.ob = glm.fit.SS(basis.x, Xt, Xv, Yt, samp.prob)

  # Supervised estimator for beta:

  beta.sl = glm.fit.ob$beta.sl

  # Our semi-supervised estimator for beta:
  beta.ssl = glm.fit.ob$beta.ssl
  beta.ssl.2 = glm.fit.ob$beta.ssl.2

  beta.true.all <- rbind(beta.true.all, beta.true)
  beta.sl.all <- rbind(beta.sl.all, beta.sl)
  beta.ssl.all <- rbind(beta.ssl.all, beta.ssl)
  beta.ssl.2.all <- rbind(beta.ssl.2.all, beta.ssl.2)
}

setwd('~/Desktop')
write.table(beta.true.all, paste0(n.t, "-beta_true.csv"),
           col.names=FALSE, quote=FALSE, row.names=FALSE)
write.table(beta.sl.all, paste0(n.t, "-beta_sl.csv"),
           col.names=FALSE, quote=FALSE, row.names=FALSE)
write.table(beta.ssl.all, paste0(n.t, "-beta_ssl.csv"),
           col.names=FALSE, quote=FALSE, row.names=FALSE)
write.table(beta.ssl.2.all, paste0(n.t, "-beta_ssl2.csv"),
           col.names=FALSE, quote=FALSE, row.names=FALSE)




