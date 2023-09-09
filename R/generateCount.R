#' Generate Mean Structure for Counts
#'
#' This function produces the mean structure based on provided parameters.
#' It is designed to be used as an internal helper function for `generateReads`.
#'
#' @param strainSizes A numeric vector specifying the size of each strain.
#' @param alpha A numeric value representing the intercept term.
#' @param beta A numeric value for the regression coefficient.
#' @param sigma A numeric matrix or value for the random effects covariance structure.
#' @param X A numeric vector or matrix of covariates.
#'
#' @return A list containing:
#'   * `means`: The computed means based on the parameters.
#'   * `num.x`: The sample size.
#'   * `col.names`: Column names formatted as "S{strain number}_{sample number within strain}".
#'
#' @keywords internal
generateMeans <- function(strainSizes, alpha, beta, sigma, X) {

  num.x <- length(X) #sample size
  num.b <- length(strainSizes) #number of strains

  # Check if sigma is univariate or multivariate and generate random effects
  if(is.null(dim(sigma))) {
    b <- stats::rnorm(num.b, sd = sqrt(sigma))
    b0 <- b
    b1 <- rep(0, num.b)
  } else {
    b <- MASS::mvrnorm(n = num.b, mu = c(0, 0), Sigma = sigma)
    b0 <- b[, 1]
    b1 <- b[, 2]
  }

  b0 <- rep(b0, strainSizes)
  b1 <- rep(b1, strainSizes)

  # compute E(y) = g^-1(alpha + beta*X + b0 + b1*X)
  means <- exp(alpha + beta * X + b0 + b1 * X)

  # name : format is S[strain number]_[sample number within strain]
  # e.g S3_2 is the second sample within strain 3
  col.names <- lapply(1:num.b, function(x) {
    paste0("S", x, "_", 1:(strainSizes[x]))
  })

  return(list(means = means, num.x = num.x, col.names = unlist(col.names)))
}

#' Generate Counts from Given Parameters and Model
#'
#' This function produces counts based on provided parameters, a mean function, and a count function.
#' It uses the `generateMeans` function to determine the mean structure and then applies a count-generating function.
#'
#' @param strainSizes A numeric vector specifying the size of each strain.
#' @param alpha A numeric value representing the intercept term.
#' @param beta A numeric value for the regression coefficient.
#' @param sigma A numeric matrix or value for the random effects covariance structure.
#' @param X A numeric vector or matrix of covariates.
#' @param phi A numeric dispersion parameter passed to the `countFun`.
#' @param countFun A function for generating count data. This function should accept parameters `n`, `mu`, and `phi` at minimum.
#' @param ... Additional arguments passed to the `countFun`.
#'
#' @return A matrix with counts. Each column corresponds to a different strain-sample and is named accordingly.
#'
#' @keywords internal
generateReads <- function(strainSizes, alpha, beta, sigma, phi, X, countFun, ...) {
  mu <- generateMeans(strainSizes, alpha, beta, sigma, X)
  #y ~ countFun(mu, phi, ...)
  counts <- lapply(1:mu$num.x, function(x) {
    countFun(n = 1, mu = mu$means[x], phi, ...)
  })
  counts <- matrix(do.call(c, counts), nrow = 1)
  colnames(counts) <- mu$col.names
  return(counts)

}

#' Generate Negative Binomial Reads
#'
#' This function produces a vector of read counts generated from a Negative Binomial distribution
#' based on the provided parameters.
#'
#' @param strainSizes A numeric vector specifying the size of each strain.
#' @param alpha A numeric specifying the intercept.
#' @param beta A numeric specifying the regression coefficient.
#' @param sigma A numeric or matrix specifying the variance (or variance-covariance matrix) for random effects.
#' @param phi A numeric specifying the dispersion parameter for the Negative Binomial distribution.
#' @param X A numeric vector or matrix of covariates.
#'
#' @return A numeric vector of Negative Binomial read counts.
#'
#' @seealso \code{\link[MASS]{rnegbin}}
#' @importFrom MASS rnegbin
#' @keywords internal
getNBReads <- function(strainSizes, alpha, beta, sigma, phi, X) {
  generateReads(strainSizes, alpha, beta, sigma, phi, X, countFun = MASS::rnegbin)
}

#' Generate Compound Poisson (Tweedie) Reads
#'
#' This function produces a vector of read counts generated from a Compound Poisson (Tweedie) distribution
#' based on the provided parameters.
#'
#' @param strainSizes A numeric vector specifying the size of each strain.
#' @param alpha A numeric specifying the intercept.
#' @param beta A numeric specifying the regression coefficient.
#' @param sigma A numeric or matrix specifying the variance (or variance-covariance matrix) for random effects.
#' @param p A numeric specifying the power parameter for the Tweedie distribution.
#' @param phi A numeric specifying the dispersion parameter for the Tweedie distribution.
#' @param X A numeric vector or matrix of covariates.
#'
#' @return A numeric vector of Compound Poisson (Tweedie) read counts.
#'
#' @seealso \code{\link[tweedie]{rtweedie}}
#' @importFrom tweedie rtweedie
#' @keywords internal
getCPReads <- function(strainSizes, alpha, beta, sigma, p, phi, X) {
  generateReads(strainSizes, alpha, beta, sigma, phi, X, countFun = tweedie::rtweedie, xi = p)
}

#' Generate a Matrix of Read Counts for Multiple Features
#'
#' This function produces a matrix of read counts for multiple features
#' based on the provided parameters and method specified (`NB` or `CP`).
#' It is designed to be an internal helper function.
#'
#' @param strainSizes A numeric vector specifying the size of each strain.
#' @param alphas A numeric vector of intercepts.
#' @param betas A numeric vector of regression coefficients.
#' @param sigma2s A list of variance-covariance matrices for random effects, one for each feature.
#' @param phis A numeric vector of dispersion parameters for the count distribution.
#' @param X A numeric vector or matrix of covariates.
#' @param method A character string. It can be either "NB" (Negative Binomial) or "CP" (some other method). Default is "NB".
#' @param ps A numeric vector. Only used if method is "CP". Default is NULL.
#'
#' @return A matrix of read counts, where each row corresponds to a feature and each column
#'   corresponds to a sample/strain.
#'
#' @keywords internal
generateReadMatrix <- function(strainSizes, alphas, betas, sigma2s, phis, X,
                               method = c("NB", "CP"), ps = NULL) {
  cat("Generating with", method)

  method <- match.arg(method)

  func <- switch(method,
                 NB = getNBReads,
                 CP = getCPReads)

  getFeatureData <- function(i) {
    if (method == "NB") {
      return(func(strainSizes, alphas[i], betas[i], sigma2s[[i]], phis[i], X))
    } else {
      return(func(strainSizes, alphas[i], betas[i], sigma2s[[i]], ps[i], phis[i], X))
    }
  }

  num.features <- length(alphas)
  countMatrix <- do.call(rbind, lapply(1:num.features, getFeatureData))
  rownames(countMatrix) <- paste("Gene ", 1:num.features)
  return(countMatrix)
}

#' Generate Read Matrix Using Negative Binomial Distribution
#'
#' This function produces a matrix of read counts generated from a Negative Binomial distribution
#' based on the provided parameters.
#'
#' @param strainSizes A numeric vector specifying the size of each strain.
#' @param alphas A numeric vector of intercepts.
#' @param betas A numeric vector of regression coefficients.
#' @param sigma2s A list of numeric values or matrices specifying the variances (or variance-covariance matrices)
#' for random effects for each gene.
#' @param phis A numeric vector specifying the dispersion parameters for the Negative Binomial distribution.
#' @param X A numeric vector or matrix of covariates.
#'
#' @return A numeric matrix of Negative Binomial read counts where each row represents a gene
#' and each column corresponds to a sample.
#'
#' @seealso \code{\link[MASS]{rnegbin}}
#' @importFrom MASS rnegbin
#' @export
getNBReadMatrix <- function(strainSizes, alphas, betas, sigma2s, phis, X)
  return(generateReadMatrix(strainSizes, alphas, betas, sigma2s, phis, X, "NB"))

#' Generate Read Matrix Using Compound Poisson Distribution
#'
#' This function produces a matrix of read counts generated from a Compound Poisson distribution
#' based on the provided parameters.
#'
#' @param strainSizes A numeric vector specifying the size of each strain.
#' @param alphas A numeric vector of intercepts.
#' @param betas A numeric vector of regression coefficients.
#' @param sigma2s A list of numeric values or matrices specifying the variances (or variance-covariance matrices)
#' for random effects for each gene.
#' @param phis A numeric vector specifying the dispersion parameters for the Compound Poisson distribution.
#' @param X A numeric vector or matrix of covariates.
#' @param ps A numeric vector specifying the power parameter for the Tweedie distribution.
#'
#' @return A numeric matrix of Compound Poisson read counts where each row represents a gene
#' and each column corresponds to a sample.
#'
#' @seealso \code{\link[tweedie]{rtweedie}}
#' @importFrom tweedie rtweedie
#' @export
getCPReadMatrix <- function(strainSizes, alphas, betas, sigma2s, phis, X, ps)
  return(generateReadMatrix(strainSizes, alphas, betas, sigma2s, phis, X, "CP", ps))

