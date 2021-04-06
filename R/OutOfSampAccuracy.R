########################################################
######### Evaluate out-of-sampe AUC, BS, OMR ###########
#### Used for the comparison in Section S5  #############
########################################################

logit_pred <- function(X.test, beta.fit){
  pred_mean <- g.logit(cbind(1, X.test) %*% beta.fit)
  pred_mean[which(is.na(pred_mean))] <- ifelse(cbind(1, X.test) %*% beta.fit > 0, 1, 0)
  return(list(class = ifelse(pred_mean > 0.5, 1, 0), pred = pred_mean))
}

Eva_class <- function(Y.test, results, weight = rep(1 / length(Y.test),
                                                    length(Y.test))){

  auc <- auc(Y.test, results$class)
  mse <- sum((Y.test - results$pred)^2 * weight)
  ae <- sum(abs(Y.test - results$class) * weight)

  return(list(auc = auc, mse = mse, ae = ae))
}
