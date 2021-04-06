library(MASS)

#' Generate data from a multivariate normal
#'
#' @param n Size of data set.
#' @param mu Mean vector.
#' @param Sigma Covariance matrix.
#' @export
#' @return A matrix of multivariate normal data.
CovariateGen <- function(n, mu, sigma){

  return(mvrnorm(n, mu, sigma))
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



#' Generate Stratification Variable
#'
#' @param covariates Covariate matrix.
#' @param num_stratum Number of stratum.
#' @export
#' @return A stratification variable
StratificationVar <- function(covariates, num_stratum = 2){
  strat_var_1 <- I(covariates[,1] + rnorm(N) < 0.5)
  strat_var_3 <- I(covariates[,3] + rnorm(N) < 0.5)

  if (num_strata == 2){
    strat_var <- ifelse(strat_var_1, 1, 0)
  }

  if (num_strata == 4){
    strat_var <- rep(0, N)
    strat_var[intersect(which(strat_var_1 == 1), which(strat_var_3 == 1))] <- 1
    strat_var[intersect(which(strat_var_1 == 1), which(strat_var_3 == 0))] <- 2
    strat_var[intersect(which(strat_var_1 == 0), which(strat_var_3 == 1))] <- 3
  }

  return(strat_var)
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
DataGeneration <- function(n_lab, n_unlab, p, rho, signal = c(1, 1, 0.5, 0.5),
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
  signal = c(signal, rep(0, p - length(signal)))

  # Covariance for covariates.
  sigma = 3*ARoneCovMat(p = p, rho = rho)

  # Covariates.
  covariates = CovariateGen(N, mu = rep(0, p), sigma)

  # Linear predictor.
  lin_pred = c(covariates %*% signal)

  # Outcome and stratification variable.
  if (model_specification == 'outcome_correct_imputation_correct'){

    strat_var <- StratificationVar(covariates, num_stratum = num_stratum)
    Y <- I(lin_pred + rlogis(N) > 0)

  }

  if (model_specification == 'outcome_incorrect_imputation_correct'){

    strat_var <- StratificationVar(covariates, num_stratum = num_stratum)

    if (num_strata == 2){
      gamma.coef <- c(signal, c(0.5, 0, 0, 0.5),
                      rep(0, 8), - 0.5,rep(0, 4), -0.5)
      basis <- ns.basis(covariates, S, 3, basis.type = 'interact')
    }
    if (num_strata == 4){
      gamma.coef <- c(signal, c(0.5, 0, 0, 0.5), rep(0, 8),
                      -0.5, rep(0, 4), -0.5, 0, 0)
      basis <- ns.basis(covariates, S, 3, basis.type = 'interact')
    }

    Y = I(c(basis %*% gamma.coef) + rlogis(N) > C)
  }


  if (model_specification == 'outcome_incorrect_imputation_incorrect'){

    strat_var <- StratificationVar(covariates, num_stratum = num_stratum)
    incorrect_signal <- c(-2, rep(0,3), -3, -3, 0, 0, rep(0,2))
    Y0 <- lin_pred + covariates[,1]^2 + covariates[,3]^2 +
      rgumbel(N, -2, 0.3)*exp(covariates%*% incorrect_signal)
    Y <- I(Y0 > 0)
  }

  if (model_specification == 'outcome_incorrect_imputation_correct_supp'){

    strat_var_1 <- I(covariates[,1] + covariates[,2] + rnorm(N) > 1.5)
    strat_var <- ifelse(strat_var_1, 1, 0)

    gamma.coef <- c(signal, c(0.5, 0, 0, 0.5), rep(0, 8),
                    -0.5, rep(0, 4), -0.5)
    basis <- ns.basis(covariates, strat_var, 3, basis.type = 'interact')

    Y0 = cbind(1, basis) %*% gamma.coef * strat_var +
      (0.8 * cbind(1, basis) %*% gamma.coef - 5) * (1 - strat_var) + rlogis(N)
    Y = I(Y0 > 1)

  }

  if (model_specification == 'outcome_incorrect_imputation_incorrect_supp'){
    B = 1.5
    C = 5

    strat_var_1 <- I(covariates[,1] + covariates[,2] + rnorm(N) > 1.5)
    strat_var <- ifelse(strat_var_1, 1, 0)

    incorrect_signal = c(-2, rep(0,3), -3, -3, 0, 0, rep(0,2))
    Y0 = (lin_pred + covariates[,1]^2 + covariates[,3]^2) * strat_var +
      (0.8 * lin_pred - 5)* (1 - strat_var) +
      rgumbel(N, -2, 0.3)*exp(c(covariates%*% incorrect_signal));
    Y = I(Y0 > 1)
  }


  if (model_specification == 'gaussian_mixture'){

    strat_var <- StratificationVar(covariates, num_stratum = num_stratum)

    mu_diff <- c(0.2, -0.2, 0.2, -0.2, 0.2, -0.2, 0.1, -0.1, rep(0, p - 8))
    sigma_diff <- autocorr.mat(p, rho = 0.3) - diag(rep(0.4, p)) +
      matrix(0.2, p, p)



    Y = rbinom(N, 1, 0.5)
    covariates = matrix(0, N, p)

    mu1 = rep(0, p)
    sigma1 = autocorr.mat(p = p, rho = rho)
    covariates[which(Y == 0), ] = CovariateGen(length(which(Y == 0)), mu1, sigma1)

    mu2 = mu1 + mu_diff
    sigma2 = Sigma1 + Sigma_diff
    covariates[which(Y == 1), ] = CovariateGen(length(which(Y == 1)), mu2, sigma2)


    covariates[which(Y == 1), 3] <- covariates[which(Y == 1), 3] + 0.12 * covariates[which(Y == 1), 3]^3
    covariates[which(Y == 1), 4] <- covariates[which(Y == 1), 4] + 0.12 * covariates[which(Y == 1), 4]^3

    covariates[which(Y == 1), 7] <- covariates[which(Y == 1), 7] + 0.12 * covariates[which(Y == 1), 7]^3
    covariates[which(Y == 1), 8] <- covariates[which(Y == 1), 8] + 0.12 * covariates[which(Y == 1), 8]^3

  }

  # Size of strata
  n_each <- n_lab / num_strata

  ##ind_lst <- vector('list', num_strata)

  if (num_strata == 2){

    # Stratum 1
    stratum_1 = which(S == 1);
    ind.1 = sample(stratum_1, size = n_each);

    # Stratum 2
    stratum_2 = which(S == 0);
    ind.2 = sample(stratum_2, size = n_each);

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
    samp.prob = c(rep(n_each / N.1, n_each), rep(n_each / N.2, n_each));
  }

  if (num_strata == 4){

    # Stratum 1
    stratum_1 = which(S == 0);
    ind.1 = sample(stratum_1, size = n_each);

    # Stratum 2
    stratum_2 = which(S == 1);
    ind.2 = sample(stratum_2, size = n_each);

    # Stratum 3
    stratum_3 = which(S == 2);
    ind.3 = sample(stratum_3, size = n_each);

    # Stratum 4
    stratum_4 = which(S == 3);
    ind.4 = sample(stratum_4, size = n_each);

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
    samp.prob = c(rep(n_each/N.1, n_each), rep(n_each/N.2, n_each),
                  rep(n_each/N.3, n_each), rep(n_each/N.4, n_each));

  }

  return(list(Covariates = Covariates, Y = Y, S = S.all, St = S.t, Sv = S.v, Xt = Xt,
              Xv = Xv, Yt = Yt, samp.prob = samp.prob, signal = signal,
              Xr = Xr, Yr = Yr, Sr = strat_var, Xvr = Xvr))

}
