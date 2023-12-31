% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateCount.R
\name{getNBReads}
\alias{getNBReads}
\title{Generate Negative Binomial Reads}
\usage{
getNBReads(strainSizes, alpha, beta, sigma, phi, X)
}
\arguments{
\item{strainSizes}{A numeric vector specifying the size of each strain.}

\item{alpha}{A numeric specifying the intercept.}

\item{beta}{A numeric specifying the regression coefficient.}

\item{sigma}{A numeric or matrix specifying the variance (or variance-covariance matrix) for random effects.}

\item{phi}{A numeric specifying the dispersion parameter for the Negative Binomial distribution.}

\item{X}{A numeric vector or matrix of covariates.}
}
\value{
A numeric vector of Negative Binomial read counts.
}
\description{
This function produces a vector of read counts generated from a Negative Binomial distribution
based on the provided parameters.
}
\seealso{
\code{\link[MASS]{rnegbin}}
}
