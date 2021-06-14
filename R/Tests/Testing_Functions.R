run_data_gen_test <- function(old_data, new_data){

  signal_test <- all.equal(old_data$b0[-1], new_data$signal)

  cov_test <- all.equal(old_data$X0, new_data$covariates)
  lab_cov_test <- all.equal(old_data$Xt, new_data$covariates_lab)
  unlab_cov_test <- all.equal(old_data$Xv, new_data$covariates_unlab)
  rs_lab_cov_test <- all.equal(old_data$Xr, new_data$covariates_random_samp)
  rs_unlab_cov_test <- all.equal(old_data$Xvr, new_data$covariates_unlab_random_samp)

  lab_outcome_test <- all.equal(old_data$Yt, new_data$Y_lab)
  outcome_test <- all.equal(old_data$Y, new_data$Y)
  rs_lab_outcome_test <- all.equal(old_data$Yr, new_data$Y_random_samp)

  lab_strat_test <- all.equal(old_data$St, new_data$S_lab)
  unlab_strat_test <- all.equal(old_data$Sv, new_data$S_unlab)
  rs_lab_strat_test <- all.equal(old_data$Sr, new_data$S_random_samp)
  strat_test <- all.equal(old_data$S, new_data$strat_var)
  samp_prob_test <- all.equal(old_data$samp.prob, new_data$samp_prob)

  return(list(signal_test = signal_test, cov_test = cov_test,
              lab_cov_test = lab_cov_test, unlab_cov_test = unlab_cov_test,
              rs_lab_cov_test = rs_lab_cov_test,
              rs_unlab_cov_test = rs_unlab_cov_test,
              lab_outcome_test = lab_outcome_test, outcome_test =outcome_test,
              rs_lab_outcome_test = rs_lab_outcome_test,
              lab_strat_test = lab_strat_test,
              unlab_strat_test = unlab_strat_test,
              rs_lab_strat_test = rs_lab_strat_test, strat_test = strat_test,
              samp_prob_test = samp_prob_test))
}
