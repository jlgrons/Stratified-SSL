% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SupervisedApparentAccuracy.R
\name{SupervisedApparentAccuracy}
\alias{SupervisedApparentAccuracy}
\title{Apparent estimates for brier score (MSE) and misclassification rate (OMR).}
\usage{
SupervisedApparentAccuracy(
  X_labeled,
  y,
  beta_SL,
  samp_prob,
  resamp_weight = NULL,
  threshold = 0.5
)
}
\arguments{
\item{X_labeled}{Covariate matrix for labeled data set.}

\item{y}{Numeric outcome vector.}

\item{samp_prob}{Numeric vector of weights.}

\item{resamp_weight}{Numeric vector of resampling weights.}

\item{threshold}{Threshold for overall misclassification rate.}
}
\value{
Supervised MSE and OMR.
}
\description{
Apparent estimates for brier score (MSE) and misclassification rate (OMR).
}
