########################################################
################ SL with random forest ##################
####### Used for the comparison in Section S5  #########
########################################################


RF_predict <- function(Yt, Xt, St, X.test, S.test, n.tree = 100){
  reg <- Xt
  test.reg <- X.test

  if (is.null(St)){

  }else{
    for (j in 1:(length(unique(St)) - 1)) {
      reg <- cbind(reg, ifelse(St == j, 1, 0))
      test.reg <- cbind(test.reg, ifelse(S.test == j, 1, 0))
    }
  }


  rf_fit <- randomForest(reg, as.factor(Yt), ntree = n.tree)
  ##print(length(test.reg[1,]))
  rf_class <- predict(rf_fit, test.reg)
  rf_class <- as.numeric(rf_class) - 1

  rf_fit_reg <- randomForest(reg, Yt, ntree = n.tree)
  rf_pred <- predict(rf_fit_reg, test.reg)

  return(list(class = rf_class, pred = rf_pred))
}
