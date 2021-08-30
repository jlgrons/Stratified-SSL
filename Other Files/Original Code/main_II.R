source('function.R')

# Both outcome and imputation models are wrong (Section 7 (iii))

# Data generation:

# Sample size of label data
n.t <- 400

# Matrix correlation
rho <- 0.4

# Sample size of unlabel data
N <- 20000
p <- 10

# Number of stratum
strata.num <- 2

# Number of folds for CV
K.fold = 3

data <- data_gen(n.t, N, p, rho, model.spec = 'II', strata_num = strata.num)
basis.type <- 'spline'
knots = 3

# Generate the true regression parameters
data.v <- data_gen(n.t, 100000-n.t, p, rho, model.spec = 'II')
X0.v <- data.v$X0
Y.v <- data.v$Y
beta.true = glm(Y.v~X0.v, family = 'binomial')$coeff

# True Brier score (BS):
true.mse = mean((Y.v - g.logit(cbind(1,X0.v) %*% beta.true))^2)

# True overall misclassification rate (OMR)
true.ae = mean(abs(Y.v - I(g.logit(cbind(1,X0.v) %*% beta.true) > 0.5) ))


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

# Beta estimated without including the sampling probability:

beta.naive = glm.fit.ob$beta.naive
gamma = glm.fit.ob$gamma

# Density ratio estimator for beta:

beta.dr = glm.fit.ob$beta.dr
proj.dr = glm.fit.ob$proj.dr


K.fold = K.fold; rep = 20;
# Cross-validated residuals
resids.cv.ob = resids.cv(basis.x, Xt, Xv, Yt, samp.prob, K.fold)
resids.beta.sl = resids.cv.ob$resids.beta.sl
resids.beta.ssl = resids.cv.ob$resids.beta.ssl
resids.beta.dr = resids.cv.ob$resids.beta.dr
resids.gamma = resids.cv.ob$resids.gamma


# Standard error estimates for supervised and SS estimates of regression parameter
se.beta.sl = se.est(beta.sl, Yt, Xt, Xv, resids.beta.sl)$se
se.beta.ssl.ob = se.est(beta.ssl, Yt, Xt, Xv, resids.gamma)
se.beta.ssl = se.beta.ssl.ob$se
A = se.beta.ssl.ob$A

se.beta.dr = se.est.dr(beta.dr, Yt, Xt, Xv, resids.beta.dr, proj.dr)


# Minimum Variance Estimator (here the component-wise optimal estimator)
# constant for stability due to high correlation between beta_SL and beta_SSL
eps.s = (n.t*(se.beta.sl^2 + se.beta.ssl^2))/(2*n.t^0.6)
beta.ssl.w.ob = min.var.est(beta.ssl, beta.sl, resids.gamma,
                            resids.beta.sl, Xt, Xv, eps.s);

# Our ensembel estimator (of the SL and SSL) for beta:
beta.ssl.w = beta.ssl.w.ob$beta

w.beta = beta.ssl.w.ob$weight
se.beta.ssl.w = beta.ssl.w.ob$se.est

# Estimated asymptotic standard error of the estimators (used for confidence interval construction):

se.beta.sl
se.beta.ssl
se.beta.ssl.w
se.beta.dr


# Apparent estimates of the model accuracy (BS and OMMR):
apparent.ob = model.eval.ap(Yt, beta.sl, beta.ssl.w,
                            gamma, beta.dr, Xt, Xv, basis.x, samp.prob)


ap.sl.mse = apparent.ob$mse.sl
ap.sl.ae = apparent.ob$ae.sl

ap.naive.mse = apparent.ob$mse.naive
ap.naive.ae = apparent.ob$ae.naive

ap.dr.mse = apparent.ob$mse.dr
ap.dr.ae = apparent.ob$ae.dr

ap.ssl.mse = apparent.ob$mse.ssl
ap.ssl.ae = apparent.ob$ae.ssl

ap.ae = c(ap.ssl.ae, ap.sl.ae, ap.dr.ae, ap.naive.ae)
ap.mse = c(ap.ssl.mse, ap.sl.mse, ap.dr.mse, ap.naive.mse)

ap.ae
ap.mse



# Cross-validated (CV) estimates of the model accuracy (BS and OMMR):
cv.ob = model.eval.cv(basis.x, Xv, Xt, Yt, samp.prob,
                      w.beta[,1], K.fold = K.fold, rep = rep)

cv.sl.mse = cv.ob$mse.sl
cv.sl.ae = cv.ob$ae.sl

cv.naive.mse = cv.ob$mse.naive
cv.naive.ae = cv.ob$ae.naive

cv.ssl.mse = cv.ob$mse.ssl
cv.ssl.ae = cv.ob$ae.ssl

cv.dr.mse = cv.ob$mse.dr
cv.dr.ae = cv.ob$ae.dr

cv.naive.mse = cv.ob$mse.naive
cv.naive.ae = cv.ob$ae.naive

cv.ae = c(cv.ssl.ae, cv.sl.ae, cv.dr.ae, cv.naive.ae)
cv.mse = c(cv.ssl.mse, cv.sl.mse, cv.dr.mse, cv.naive.mse)

# Ensemble of apparent and CV estimator (Section 4)
# in the order of SSL (our method), SL, DR and SL method not considering the stratified sampling
w <- K.fold / (2 * K.fold - 1)
w * ap.ae + (1 - w) * cv.ae
w * ap.mse + (1 - w) * cv.mse

# Standard error estimation for the estimators of BS and OMR (used for CI construction)
b = 500;
pert.ob = model.eval.se(b, Yt, Xt, Xv, basis.x, samp.prob, w.beta,
                        beta.sl, beta.ssl.w, beta.dr, resids.beta.sl,
                        resids.gamma, resids.beta.dr, A, c = 0.5)


ssl.pert.mse = pert.ob$ssl.pert.mse
sl.pert.mse = pert.ob$sl.pert.mse
dr.pert.mse = pert.ob$dr.pert.mse

ssl.pert.ae = pert.ob$ssl.pert.ae
sl.pert.ae = pert.ob$sl.pert.ae
dr.pert.ae = pert.ob$dr.pert.ae

sd(sl.pert.mse)/sd(ssl.pert.mse)
sd(sl.pert.ae)/sd(ssl.pert.ae)

sd(sl.pert.ae)/sd(dr.pert.ae)
sd(sl.pert.mse)/sd(dr.pert.mse)

sd(sl.pert.mse)*100
sd(sl.pert.ae)*100
sd(ssl.pert.mse)*100
sd(ssl.pert.ae)*100



# Evaluate out-of-sample Classification and Prediction Performance (for results in Section S5)
# Generate test data

data.test <- data_gen(n.t = 1000, n.v = 9000, p, rho, model.spec = 'II', strata_num = strata.num)

X.test <- data.test$X0
Y.test <- data.test$Y
S.test <- data.test$S

# Random Forest

rf_model <- RF_predict(Yt, Xt, St, X.test, S.test)

# AUC, BS and OMR.
eva.rf <- Eva_class(Y.test, rf_model)
eva.sl <- Eva_class(Y.test, logit_pred(X.test, beta.sl))
eva.ssl <- Eva_class(Y.test, logit_pred(X.test, beta.ssl))
eva.dr <- Eva_class(Y.test, logit_pred(X.test, beta.dr))
eva.ssl.w <- Eva_class(Y.test, logit_pred(X.test, beta.ssl.w))

auc.perform <- c(eva.rf$auc, eva.sl$auc,
                 eva.ssl$auc, eva.dr$auc, eva.ssl.w$auc)

ae.perform <- c(eva.rf$ae, eva.sl$ae,
                eva.ssl$ae, eva.dr$ae, eva.ssl.w$ae)

mse.perform <- c(eva.rf$mse, eva.sl$mse,
                 eva.ssl$mse, eva.dr$mse, eva.ssl.w$mse)
auc.perform
ae.perform
mse.perform
