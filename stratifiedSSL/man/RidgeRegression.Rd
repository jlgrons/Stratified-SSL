% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RidgeRegression.R
\name{RidgeRegression}
\alias{RidgeRegression}
\title{Ridge regression.}
\usage{
RidgeRegression(X, y, weights = NULL, lambda = 1e-04, family = "binomial")
}
\arguments{
\item{X}{Covariate matrix.}

\item{y}{Numeric outcome vector.}

\item{weights}{Numeric vector of weights.}

\item{lambda}{Penalization parameter.}

\item{family}{Exponential family of interest.}
}
\value{
Vector containing regression coefficients.
}
\description{
Ridge regression.
}
