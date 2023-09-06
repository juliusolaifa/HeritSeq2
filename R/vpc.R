computeVPC <- function(para, covar, type) {

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
                   sig22=0, sig12=0, phi=row[length(row)], x=covar)

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

computeNBVPC <- function(para) {
  slope_check <- ifelse(is.null(dim(para)), length(para), dim(para)[2])
  if(slope_check == 4) {
    computeVPC(para,0,"NB")
  } else {
    computeVPC(para,1,"NB")
  }
}

computeCPVPC <- function(para) {
  slope_check <- ifelse(is.null(dim(para)), length(para), dim(para)[2])
  if(slope_check == 5) {
    computeVPC(para,0,"CP")
  } else {
    computeVPC(para,1,"CP")
  }
}
