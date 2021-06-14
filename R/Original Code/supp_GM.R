
# Gaussian mixture (GM) setting in Section S5 = 'GM'

source('function.R')

##### Both outcome and imputation models are correct (Section 7 (i)) #####

# Data generation:

# Sample size of label data
n.t <- 200

# Matrix correlation
rho = 0.2

# Sample size of unlabel data
N <- 20000
p <- 10

# Number of stratum
strata.num <- 2

# Number of folds for CV
K.fold = 3
basis.type = 'interact'
knots = 3


# Estimate the projection matrices of principle component analysis (PCA)

data.unsup <- data_gen(n.t, N - n.t, p, rho,
                       model.spec = 'GM', strata_num = strata.num)
out.num <- 5
basis.num <- 5
data.pca <- extract.pc(data.unsup$X0, out.num = out.num, basis.num = basis.num)
V.first = data.pca$V.first
V.second = data.pca$V.second

# Generate the true parameters
data.v <- data_gen(n.t, 100000 - n.t, p, rho, model.spec = 'GM', strata_num = strata.num)
X0.v <- data.v$X0
Y.v <- data.v$Y
data.pca.v <- extract.pc(X0.v, out.num = out.num, basis.num = basis.num,
                         V.first = V.first, V.second = V.second)
X.v <- data.pca.v$x

# True regression parameters:
beta.true = glm(Y.v ~ X.v, family = 'binomial')$coeff

# True Brier score (BS):
mse.true = mean((Y.v - g.logit(cbind(1, X.v) %*% beta.true))^2)

# True overall misclassification rate (OMR)
ae.true = mean(abs(Y.v - I(g.logit(cbind(1, X.v) %*% beta.true) > 0.5)))


# Generate the data for training and in-sample model accuracy evaluation:

data <- data_gen(n.t, N, p, rho, model.spec = 'GM', strata_num = strata.num)
Xt <- data$Xt
Yt <- data$Yt
Xv <- data$Xv
samp.prob <- data$samp.prob
Strata_ID <- data$S
St <- data$St


# Generate test data
data.test <- data_gen(n.t = 500, n.v = 49500, p, rho,
                      model.spec = 'GM', strata_num = strata.num)

X.test <- data.test$X0
X.test.raw <- X.test
Y.test <- data.test$Y
S.test <- data.test$S
basis.test.raw <- ns.basis(X.test, S.test, knots, basis.type = basis.type)


# Train random forest (for purpose of prediction):

rf_model_origin <- RF_predict(Yt, Xt, St, X.test, S.test, n.tree = 100)
eva.rf.origin <- Eva_class(Y.test, rf_model_origin)

# Performance of random forest

eva.rf.origin

# Transform the test data with the PCA matrices and form the basis function of imputation model.

basis.test <- ns.basis(X.test, S.test, knots, basis.type = basis.type)
data.pca.test <- extract.pc(X.test, out.num = out.num, basis.num = basis.num,
                            V.first = V.first, V.second = V.second)
X.test <- data.pca.test$x
basis.test <- cbind(X.test, basis.test)

weights = 1/samp.prob/mean(1/samp.prob)
dat.all = rbind(Xt, Xv)

### Basis of the raw data

basis.x.raw <- ns.basis(dat.all, Strata_ID, knots, basis.type = basis.type)

## Performance of the SL estimator with the design of raw data

beta.raw.sl = glm(Yt ~ Xt, family = 'binomial', weights = weights)$coeff;
eva.raw.sl <- Eva_class(Y.test, logit_pred(X.test.raw, beta.raw.sl))

eva.raw.sl

fit.raw.ssl <- glm.fit.SS(basis.x.raw, Xt, Xv, Yt, samp.prob)
beta.raw.ssl <- fit.raw.ssl$beta.ssl
gamma.raw <- fit.raw.ssl$gamma

## Performance of the imputation model of the basis of raw data

eva.imp.raw <- Eva_class(Y.test, logit_pred(basis.test.raw, gamma.raw))
eva.imp.raw

## Performance of the SSL estimator with the design of raw data

eva.raw.ssl <- Eva_class(Y.test, logit_pred(X.test.raw, beta.raw.ssl))
eva.raw.ssl


#### Methods fitted with the PCA transformed design:

X.all <- rbind(Xt, Xv)
basis.x = ns.basis(X.all, Strata_ID, knots, basis.type = basis.type)

data.pca <- extract.pc(X.all, out.num = out.num, basis.num = basis.num,
                       V.first = V.first, V.second = V.second)

Xt <- data.pca$x[1:n.t,]
Xv <- data.pca$x[(n.t + 1):(N + n.t),]
dat.all = rbind(Xt, Xv)
basis.x <- cbind(dat.all, basis.x[,-c(1:10)])

glm.fit.ob = glm.fit.SS(basis.x, Xt, Xv, Yt, samp.prob)
beta.sl = glm.fit.ob$beta.sl
beta.ssl = glm.fit.ob$beta.ssl
gamma = glm.fit.ob$gamma
beta.dr = glm.fit.ob$beta.dr
proj.dr = glm.fit.ob$proj.dr


K.fold = K.fold; rep = 10;
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
beta.ssl.w = beta.ssl.w.ob$beta
w.beta = beta.ssl.w.ob$weight
se.beta.ssl.w = beta.ssl.w.ob$se.est


## Performance of the estimators with the design tranformed through PCA

eva.sl <- Eva_class(Y.test, logit_pred(X.test, beta.sl))
eva.dr <- Eva_class(Y.test, logit_pred(X.test, beta.dr))
eva.ssl.w <- Eva_class(Y.test, logit_pred(X.test, beta.ssl.w))
eva.ssl <- Eva_class(Y.test, logit_pred(X.test, beta.ssl))


# Apparent estimates of the model accuracy (BS and OMR):
apparent.ob = model.eval.ap(Yt, beta.sl, beta.ssl,
                            gamma, beta.dr, Xt, Xv, basis.x, samp.prob)

ap.sl.mse = apparent.ob$mse.sl
ap.sl.ae = apparent.ob$ae.sl

ap.dr.mse = apparent.ob$mse.dr
ap.dr.ae = apparent.ob$ae.dr

ap.ssl.mse = apparent.ob$mse.ssl
ap.ssl.ae = apparent.ob$ae.ssl

ap.ae = c(ap.sl.ae, ap.ssl.ae, ap.dr.ae)
ap.mse = c(ap.sl.mse, ap.ssl.mse, ap.dr.mse)


# Cross-validated (CV) estimates of the model accuracy (BS and OMR):
cv.ob = model.eval.cv(basis.x, Xv, Xt, Yt, samp.prob, w.beta[,1],
                      K.fold = K.fold, rep = rep)

cv.sl.mse = cv.ob$mse.sl
cv.sl.ae = cv.ob$ae.sl

cv.ssl.mse = cv.ob$mse.ssl
cv.ssl.ae = cv.ob$ae.ssl

cv.dr.mse = cv.ob$mse.dr
cv.dr.ae = cv.ob$ae.dr

cv.ae = c(cv.sl.ae, cv.ssl.ae, cv.dr.ae)
cv.mse = c(cv.sl.mse, cv.ssl.mse, cv.dr.mse)

# Ensemble of apparent and CV estimator (Section 4)
# in the order of SSL (our method), SL, DR and SL method not considering the stratified sampling
w <- K.fold / (2 * K.fold - 1)
w * ap.ae + (1 - w) * cv.ae
w * ap.mse + (1 - w) * cv.mse

# Evaluating them using the true BS and OMR

mse.true
ae.true
