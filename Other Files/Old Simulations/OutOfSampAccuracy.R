########################################################
######### Evaluate out-of-sampe AUC, BS, OMR ###########
#### Used for the comparison in Section S5  ############
########################################################

pred_prob <- function(X_test, beta, threshold = 0.5){
  pred_mean <- Expit(cbind(1, X_test) %*% beta)

  # Ask Molei about this.
  pred_mean[which(is.na(pred_mean))] <- ifelse(cbind(1, X_test) %*% beta > 0, 1, 0)

  return(list(class = ifelse(pred_mean > threshold, 1, 0),
              pred = pred_mean))
}


classification_evaluation <- function(y_test, results,
                                      weight = rep(1 / length(y_test),
                                                   length(y_test))){

  auc <- auc(y_test, results$class)
  mse <- sum((y_test - results$pred)^2 * weight)
  ae <- sum(abs(y_test - results$class) * weight)

  return(list(auc = auc, mse = mse, ae = ae))
}
