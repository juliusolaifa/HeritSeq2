generateMeans <- function(reps, alpha, beta, sigma, X) {

  num.x <- length(X) #sample size
  num.b <- length(reps) #number of strains

  # generate the random effect (uni-variate or multi-variate)
  b <- if(is.null(dim(sigma))) stats::rnorm(num.b, sd = sqrt(sigma))
  else MASS::mvrnorm(n=num.b, mu=c(0,0), Sigma=sigma)

  # assign column 1 to b0 and column 2 to b1 form matrix b
  b0 <- if (is.null(dim(sigma))) b else b[,1]
  b1 <- if (is.null(dim(sigma))) rep.int(0, num.b) else b[,2]
  b0 <- rep.int(b0, reps)
  b1 <- rep.int(b1, reps)

  # compute E(y) = g^-1(alpha + beta*X + b0 + b1*X)
  means <- exp(alpha + beta * X + b0 + b1 * X)

  # name : format is S[strain number]_[sample number within strain]
  # e.g S3_2 is the second sample within strain 3
  col.names <- lapply(1:num.b, function(x) {
    paste0("S", x, "_", 1:(reps[x]))
  })
  col.names <- as.vector(do.call(c, col.names))

  return(list(means = means, num.x = num.x, col.names = col.names))
}

generateReads <- function(reps, alpha, beta, sigma, X, phi, countFun, ...) {
  mu <- generateMeans(reps, alpha, beta, sigma, X)
  # y ~ countFun(mu, phi, ...)
  counts <- lapply(1:mu$num.x, function(x) {
    countFun(n = 1, mu = mu$means[x], phi, ...)
  })
  counts <- matrix(do.call(c, counts), nrow = 1)
  colnames(counts) <- mu$col.names
  return(counts)
}

getNBReads <- function(reps, alpha, beta, sigma, phi, X) {
  generateReads(reps, alpha, beta, sigma, X, phi, countFun = MASS::rnegbin)
}

getCPReads <- function(reps, alpha, beta, sigma, p, phi, X) {
  generateReads(reps, alpha, beta, sigma, X, phi, countFun = tweedie::rtweedie, xi = p)
}

generateReadMatrix <- function(reps, alphas, betas, sigma2s, phis, X,
                               method = c("NB", "CP"), ps = NULL) {
  print(paste("Generating with", method))

  method <- match.arg(method)

  func <- switch(method,
                 NB = getNBReads,
                 CP = getCPReads)

  getFeatureData <- function(i) {
    if (method == "NB") {
      return(func(reps, alphas[i], betas[i], sigma2s[[i]], phis[i], X))
    } else {
      return(func(reps, alphas[i], betas[i], sigma2s[[i]], ps[i], phis[i], X))
    }
  }

  num.features <- length(alphas)
  countMatrix <- do.call(rbind, lapply(1:num.features, getFeatureData))
  rownames(countMatrix) <- paste("Gene ", 1:num.features)
  return(countMatrix)
}

getNBReadMatrix <- function(reps, alphas, betas, sigma2s, phis, X)
  return(generateReadMatrix(reps, alphas, betas, sigma2s, phis, X, "NB"))

getCPReadMatrix <- function(reps, alphas, betas, sigma2s, phis, X, ps)
  return(generateReadMatrix(reps, alphas, betas, sigma2s, phis, X, "CP", ps))

