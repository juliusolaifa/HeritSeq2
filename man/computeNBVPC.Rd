% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vpc.R
\name{computeNBVPC}
\alias{computeNBVPC}
\title{Compute the Variance Partition Coefficient (VPC) for Negative Binomial Model}
\usage{
computeNBVPC(para, group)
}
\arguments{
\item{para}{A matrix or vector of model parameters for the Negative Binomial model.}

\item{group}{a scalar 0, 1 for binary group or some float for continuous variable}
}
\value{
A matrix with computed VPC values for the Negative Binomial model. The number of rows matches the number of rows of the input \code{para}.
}
\description{
Computes the VPC for a given set of parameters for the Negative Binomial model.
It checks the dimension of the input parameters to determine the type of computation (with slope or without).
}
\examples{
# Sample use of computeNBVPC function
sample_params <- matrix(rnorm(20), ncol=4)
computeNBVPC(sample_params,0)

}
\seealso{
\code{\link{computeVPC}}
}
