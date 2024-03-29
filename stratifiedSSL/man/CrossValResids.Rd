% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CrossValResids.R
\name{CrossValResids}
\alias{CrossValResids}
\title{Apparent estimates for brier score (MSE) and misclassification rate (OMR).}
\usage{
CrossValResids(
  basis_labeled,
  basis_unlabeled,
  X_labeled,
  X_unlabeled,
  y,
  samp_prob,
  num_folds,
  lambda = NULL
)
}
\arguments{
\item{basis_labeled}{Basis matrix for labeled data set.}

\item{X_labeled}{Covariate matrix for labeled data set.}

\item{X_unlabeled}{Covariate matrix for unlabeled data set.}

\item{y}{Numeric outcome vector.}

\item{samp_prob}{Numeric vector of weights.}

\item{num_folds}{Number of folds for cross-validation.}

\item{lambda}{Regularization parameter for the imputation model.}

\item{bais_unlabeled}{Basis matrix for unlabeled data set.}
}
\value{
Cross-validated residuals.
}
\description{
Apparent estimates for brier score (MSE) and misclassification rate (OMR).
}
