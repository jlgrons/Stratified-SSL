########################
### Helper Functions ###
########################

# Necessary packages
library(matrixStats)
library(stepPlr)
library(evd)
library(methods)
library(MASS)
library(glmnet)
library(randomForest)



# Useful functions for logistic model
dg.logit = function(xx){
  ddd <- exp(xx)/(exp(xx)+1)^2
  ddd[which(is.na(ddd))] = 0
  return(ddd)
}

logit = function(xx){log((xx)/(1-xx))};

g.logit = function(xx){1/(exp(- xx) + 1)}
