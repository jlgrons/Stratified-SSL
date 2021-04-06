library(MASS)

#' Generate data from a multivariate normal
#'
#' @param n Size of data set.
#' @param mu Mean vector.
#' @param Sigma Covariance matrix.
#' @export
#' @return A matrix of multivariate normal data.
CovariateGen <- function(n, mu, Sigma){

  return(mvrnorm(n, mu, Sigma))
}


#' Generate AR-1 Covariance Matrix
#'
#' @param p Dimension of matrix.
#' @param rho Covariance parameter.
#' @export
#' @return An AR-1 covariance matrix.
ARoneCovMat <- function(p, rho){

  id_matrix <- diag(p)
  ar_one_matrix <- rho^abs(row(id_matrix)-col(id_matrix))

  return(ar_one_matrix)
}


#' Generate AR-1 Covariance Matrix
#'
#' @param n_lab Size of labeled data.
#' @param n_unlab Size of unlabeled data.
#' @param p Number of covariates.
#' @param rho Covariance parameter for AR-1 covariance structure for covariates.
#' @param signal Values of nonzero regression parameters.
#' @param model_specification Choice of model specification.
#' @param num_strat Number of strata for stratified sampling.
#' @export
#' @return List of relevant data objects.
DataGeneration <- function(n_lab, n_unlab, p, rho, signal = c(1, 1, 0.5, 0.,5),
                     model_specification = 'CC', num_strata = 2){

  ## model_specification:
  ## 'outcome_correct_imp_correct'  = 'CC'
  ## 'outcome_wrong_imp_correct'  = 'IC'
  ## 'outcome_wrong_imp_wrong'  = 'II'
  ## Supplement:
  ## 'outcome_wrong_imp_correct'  = 'IC1'
  ## 'outcome_wrong_imp_wrong'  = 'II1'
  ## Gaussian mixture (GM) setting in Section S5 = 'GM')

  # Total data size.
  N = n_lab + n_unlab

  # Regression parameter.
  signal = c(0, signal, rep(0, p - length(signal)));

  # Covariance for covariates.
  Sigma = 3*ARoneCovMat(p = p, rho = rho)

  # Covariates.
  Covariates = CovariateGen(N, mu = rep(0, p), Sigma);

  # Linear predictor
  lin.pred = c(cbind(1,Covariates) %*% signal);

  # Outcome - Correct 'C' or Mis-specified Model otherwise
  if (model_specification == 'CC'){
    B = 0.5
    C = 0

    S1 <- I(Covariates[,1] + rnorm(N, 0, 1) < B)
    S3 <- I(Covariates[,3] + rnorm(N, 0, 1) < B)

    if (num_strata == 2){
      S <- ifelse(S1, 1, 0)
    }

    if (num_strata == 4){
      S <- rep(0, N)
      S[intersect(which(S1 == 1), which(S3 == 1))] <- 1
      S[intersect(which(S1 == 1), which(S3 == 0))] <- 2
      S[intersect(which(S1 == 0), which(S3 == 1))] <- 3
    }
    Y0 = lin.pred + rlogis(N)
    Y = I(Y0 > C)
  }

  if (model_specification == 'IC'){

    B = 0.5
    C = 0

    S1 <- I(Covariates[,1] + rnorm(N, 0, 1) < B)
    S3 <- I(Covariates[,3] + rnorm(N, 0, 1) < B)

    if (num_strata == 2){
      S <- ifelse(S1, 1, 0)
    }

    if (num_strata == 4){
      S <- rep(0, N)
      S[intersect(which(S1 == 1), which(S3 == 1))] <- 1
      S[intersect(which(S1 == 1), which(S3 == 0))] <- 2
      S[intersect(which(S1 == 0), which(S3 == 1))] <- 3
    }

    if (num_strata == 2){
      gamma.coef <- c(signal, c(0.5, 0, 0, 0.5), rep(0, 8), - 0.5, rep(0, 4), - 0.5)
      basis <- ns.basis(Covariates, S, 3, basis.type = 'interact')
    }
    if (num_strata == 4){
      gamma.coef <- c(signal, c(0.5, 0, 0, 0.5), rep(0, 8), - 0.5, rep(0, 4), - 0.5, 0, 0)
      basis <- ns.basis(Covariates, S, 3, basis.type = 'interact')
    }

    Y0 = cbind(1, basis) %*% gamma.coef + rlogis(N);
    Y = I(Y0 > C)
  }


  if (model_specification == 'II'){
    B = 0.5
    C = 0

    S1 <- I(Covariates[,1] + rnorm(N, 0, 1) < B)
    S3 <- I(Covariates[,3] + rnorm(N, 0, 1) < B)

    if (num_strata == 2){
      S <- ifelse(S1, 1, 0)
    }

    if (num_strata == 4){
      S <- rep(0, N)
      S[intersect(which(S1 == 1), which(S3 == 1))] <- 1
      S[intersect(which(S1 == 1), which(S3 == 0))] <- 2
      S[intersect(which(S1 == 0), which(S3 == 1))] <- 3
    }
    c0 = c(0, -2, rep(0,3), -3, -3, 0, 0, rep(0,2))
    Y0 = lin.pred + Covariates[,1]^2 + Covariates[,3]^2 + rgumbel(N, -2, 0.3)*exp(cbind(1,Covariates)%*% c0);
    Y = I(Y0 > C)
  }

  if (model_specification == 'IC1'){
    B = 1.5
    C = 5

    S1 <- I(Covariates[,1] + Covariates[,2] + rnorm(N, 0, 1) > B)
    if (num_strata == 2){
      S <- ifelse(S1, 1, 0)
    }

    gamma.coef <- c(signal, c(0.5, 0, 0, 0.5), rep(0, 8), - 0.5, rep(0, 4), - 0.5)
    basis <- ns.basis(Covariates, S, 3, basis.type = 'interact')

    Y0 = cbind(1, basis) %*% gamma.coef * S +
      (0.8 * cbind(1, basis) %*% gamma.coef - C) * (1 - S) + rlogis(N);
    Y = I(Y0 > 1)
  }

  if (model_specification == 'II1'){
    B = 1.5
    C = 5

    S1 <- I(Covariates[,1] + Covariates[,2] + rnorm(N, 0, 1) > B)
    if (num_strata == 2){
      S <- ifelse(S1, 1, 0)
    }

    c0 = c(0, -2, rep(0,3), -3, -3, 0, 0, rep(0,2))
    Y0 = (lin.pred + Covariates[,1]^2 + Covariates[,3]^2) * S + (0.8 * lin.pred - C)* (1 - S) +
      rgumbel(N, -2, 0.3)*exp(cbind(1,Covariates)%*% c0);
    Y = I(Y0 > 1)
  }


  if (model_specification == 'GM'){
    B = 0.5
    mu_diff <- c(0.2, -0.2, 0.2, -0.2, 0.2, -0.2, 0.1, -0.1, rep(0, p - 8))
    Sigma_diff <- autocorr.mat(p, rho = 0.3) - diag(rep(0.4, p)) + matrix(0.2, p, p)

    S1 <- I(Covariates[,1] + rnorm(N, 0, 1) < B)
    S3 <- I(Covariates[,3] + rnorm(N, 0, 1) < B)

    if (num_strata == 2){
      S <- ifelse(S1, 1, 0)
    }

    if (num_strata == 4){
      S <- rep(0, N)
      S[intersect(which(S1 == 1), which(S3 == 1))] <- 1
      S[intersect(which(S1 == 1), which(S3 == 0))] <- 2
      S[intersect(which(S1 == 0), which(S3 == 1))] <- 3
    }

    Y = rbinom(N, 1, 0.5)
    Covariates = matrix(0, N, p)

    mu1 = rep(0, p)
    Sigma1 = autocorr.mat(p = p, rho = rho)
    Covariates[which(Y == 0), ] = CovariateGen(length(which(Y == 0)), mu1, Sigma1)

    mu2 = mu1 + mu_diff
    Sigma2 = Sigma1 + Sigma_diff
    Covariates[which(Y == 1), ] = CovariateGen(length(which(Y == 1)), mu2, Sigma2)

    S1 <- I(Covariates[,1] + rnorm(N, 0, 1) < B)
    S3 <- I(Covariates[,3] + Covariates[,4] + rnorm(N, 0, 1) > B)

    if (num_strata == 2){
      S <- ifelse(S1, 1, 0)
    }

    if (num_strata == 4){
      S <- rep(0, N)
      S[intersect(which(S1 == 1), which(S3 == 1))] <- 1
      S[intersect(which(S1 == 1), which(S3 == 0))] <- 2
      S[intersect(which(S1 == 0), which(S3 == 1))] <- 3
    }
    Covariates[which(Y == 1), 3] <- Covariates[which(Y == 1), 3] + 0.12 * Covariates[which(Y == 1), 3]^3
    Covariates[which(Y == 1), 4] <- Covariates[which(Y == 1), 4] + 0.12 * Covariates[which(Y == 1), 4]^3

    Covariates[which(Y == 1), 7] <- Covariates[which(Y == 1), 7] + 0.12 * Covariates[which(Y == 1), 7]^3
    Covariates[which(Y == 1), 8] <- Covariates[which(Y == 1), 8] + 0.12 * Covariates[which(Y == 1), 8]^3

  }

  # Size of strata

  n.each <- n_lab / num_strata

  ##ind_lst <- vector('list', num_strata)

  if (num_strata == 2){

    # Stratum 1
    stratum.1 = which(S == 1);
    ind.1 = sample(stratum.1, size = n.each);

    # Stratum 2
    stratum.2 = which(S == 0);
    ind.2 = sample(stratum.2, size = n.each);

    ### Labeled L & Unlabeled U datasets ###

    # Indices of labeled and unlabled data
    ind.t = c(ind.1, ind.2);
    ind.v = setdiff(1:N, ind.t);

    # Full data
    dat.N = cbind(Y, Covariates);

    # Labeled and Unlabed data sets

    ## Stratified sampling

    dat.t = dat.N[ind.t, ];
    S.t = S[ind.t]
    dat.v = dat.N[ind.v, ];
    S.v = S[ind.v]
    S.all = c(S.t, S.v)

    ## Random sampling (only used under our settings in Section S4)

    dat.r = dat.N[1:n_lab, ]
    Yr = dat.r[,1];
    Xr = dat.r[,-1];
    Xvr = dat.N[-c(1:n_lab), -1];

    # Relevant labeled data
    Yt = dat.t[,1];
    Xt = dat.t[,-1];

    # Relevant unlabeled data
    Xv = dat.v[,-1];

    # Sampling probabilities
    ind.S1.v = (which(S > 0));
    ind.S2.v  = (which(S <= 0));
    N.1 = length(ind.S1.v)
    N.2 = length(ind.S2.v)
    samp.prob = c(rep(n.each / N.1, n.each), rep(n.each / N.2, n.each));
  }

  if (num_strata == 4){

    # Stratum 1
    stratum.1 = which(S == 0);
    ind.1 = sample(stratum.1, size = n.each);

    # Stratum 2
    stratum.2 = which(S == 1);
    ind.2 = sample(stratum.2, size = n.each);

    # Stratum 3
    stratum.3 = which(S == 2);
    ind.3 = sample(stratum.3, size = n.each);

    # Stratum 4
    stratum.4 = which(S == 3);
    ind.4 = sample(stratum.4, size = n.each);

    ### Labeled L & Unlabeled U datasets ###

    # Indices of labeled and unlabled data
    ind.t = c(ind.1, ind.2, ind.3, ind.4);
    ind.v = setdiff(1:N, ind.t);

    # Full data
    dat.N = cbind(Y, Covariates);

    # Labeled and Unlabed data sets
    dat.t = dat.N[ind.t, ];
    S.t = S[ind.t]
    dat.v = dat.N[ind.v, ];
    S.v = S[ind.v]
    S.all = c(S.t, S.v)

    ## Random sampling (only used under our settings in Section S4)

    dat.r = dat.N[1:n_lab, ]
    Yr = dat.r[,1];
    Xr = dat.r[,-1];
    Xvr = dat.N[-c(1:n_lab), -1];

    # Relevant labeled data
    Yt = dat.t[,1];
    Xt = dat.t[,-1];

    # Relevant unlabeled data
    Xv = dat.v[,-1];

    # Sampling probabilities
    ind.S1.v = (which(S == 0));
    ind.S2.v  = (which(S == 1));
    ind.S3.v = (which(S == 2));
    ind.S4.v  = (which(S == 3));
    N.1 = length(ind.S1.v)
    N.2 = length(ind.S2.v)
    N.3 = length(ind.S3.v)
    N.4 = length(ind.S4.v)
    samp.prob = c(rep(n.each/N.1, n.each), rep(n.each/N.2, n.each),
                  rep(n.each/N.3, n.each), rep(n.each/N.4, n.each));

  }

  return(list(Covariates = Covariates, Y = Y, S = S.all, St = S.t, Sv = S.v, Xt = Xt,
              Xv = Xv, Yt = Yt, samp.prob = samp.prob, signal = signal,
              Xr = Xr, Yr = Yr, Sr = S, Xvr = Xvr))

}
