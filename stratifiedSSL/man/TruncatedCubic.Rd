% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelperFunctions.R
\name{TruncatedCubic}
\alias{TruncatedCubic}
\title{Truncated Cubic Function}
\usage{
TruncatedCubic(x, knot_location)
}
\arguments{
\item{x}{Numeric vector of interest.}

\item{knot_location}{Knot location for truncation.}
}
\value{
Numeric vector.
}
\description{
Calculates \eqn{((X > x)(X-x))^3}.
}
