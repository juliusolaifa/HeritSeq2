#' Compute the Variance Partition Coefficient (VPC)
#'
#' An internal function that computes the VPC for given parameters and a covariateiate.
#' This function can handle both single parameter vectors and matrices.
#'
#' @param para A matrix or vector of model parameters.
#' @param covariate A numeric vector of covariateiate values corresponding to the rows of `para`.
#' @param type Character string specifying the type of model ('NB' for Negative Binomial or 'CP' for Compound Poisson).
#'
#' @return A matrix with computed VPC values. The number of rows matches the number of rows of the input `para`.
#'
#' @keywords internal
computeVPC <- function(para, covariate, type) {

  # Inner function to compute VPC for a single row of parameters
  compute1VPC <- function(alpha, beta, sig11, sig22, sig12, phi, x, power = NULL) {
    if(is.na(sig11) || sig11 < 0 || sig22 < 0 || phi <= 0) {
      return(NA)
    }

    temp <- sig11 + 2*sig12*x + sig22 * x^2
    exp_temp <- exp(temp)

    den <- if(is.null(power)) {
      exp_temp - 1 + phi*exp_temp + exp(-alpha-beta*x-temp/2)
    } else {
      exp_temp - 1 + phi*exp((power-2)*(alpha+beta*x) + ((power^2/2)-1)*temp)
    }

    if(den == 0) {
      return(NA)
    }

    return((exp_temp - 1) / den)
  }

  # Helper function to compute VPC for each row in the matrix
  compute_vpc <- function(row) {
    params <- list(alpha=row[1], beta=row[2], sig11=row[3],
                   sig22=0, sig12=0, phi=row[length(row)], x=covariate)

    num_params <- length(row)

    # Update additional parameters based on type and number of parameters
    if (type == "CP") {
      params$power <- row[length(row)-1]
    }
    if(num_params > 4) {
      params$sig22 <- row[4]
      params$sig12 <- row[5]
    }

    return(do.call(compute1VPC, params))
  }

  # Compute VPCs for the entire dataset (either a single row or a matrix)
  if(is.null(dim(para))) {
    vpcs <- compute_vpc(para)
  } else {
    vpcs <- apply(para, 1, compute_vpc)
  }

  # Format and return the results
  vpcs <- matrix(vpcs, ncol = 1)
  rownames(vpcs) <- rownames(para)
  return(vpcs)
}

#' Compute the Variance Partition Coefficient (VPC) for Negative Binomial Model
#'
#' Computes the VPC for a given set of parameters for the Negative Binomial model.
#' It checks the dimension of the input parameters to determine the type of computation (with slope or without).
#'
#' @param para A matrix or vector of model parameters for the Negative Binomial model.
#' @param group a scalar 0, 1 for binary group or some float for continuous variable
#' @return A matrix with computed VPC values for the Negative Binomial model. The number of rows matches the number of rows of the input `para`.
#' @export
#' @examples
#' # Sample use of computeNBVPC function
#' sample_params <- matrix(rnorm(20), ncol=4)
#' computeNBVPC(sample_params,0)
#'
#' @seealso \code{\link{computeVPC}}
computeNBVPC <- function(para, group=0) {
  print(paste("Estimating VPC with NB for group", group))
  computeVPC(para,group,"NB")
}

#' Compute the Variance Partition Coefficient (VPC) for Compound Poisson Model
#'
#' Computes the VPC for a given set of parameters for the Compound Poisson model.
#' It checks the dimension of the input parameters to determine the type of computation (with slope or without).
#'
#' @param para A matrix or vector of model parameters for the Compound Poisson model.
#' @param group a scalar 0, 1 for binary group or some float for continuous variable
#' @return A matrix with computed VPC values for the Compound Poisson model. The number of rows matches the number of rows of the input `para`.
#' @export
#' @examples
#' # Sample use of computeCPVPC function
#' sample_params <- matrix(rnorm(25), ncol=5)
#' computeCPVPC(sample_params,0)
#'
#' @seealso \code{\link{computeVPC}}
computeCPVPC <- function(para,group=0) {
  cat("Estimating VPC with CP for group", group)
  computeVPC(para,group,"CP")
}
