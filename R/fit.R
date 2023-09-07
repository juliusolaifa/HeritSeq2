#' Extract Model Parameters
#'
#' Internal function to extract fixed effects, variance components, and other parameters
#' from a glmmTMB model.
#'
#' @param model A fitted \code{glmmTMB} model object.
#' @param type A character string specifying the type of the model, either "NB" for negative binomial
#' or "CP" for compound Poisson.
#' @param slope A logical indicating whether the model includes a slope for the random effect.
#' Default is FALSE.
#'
#' @return A numeric vector of extracted model parameters.
#'
#' @keywords internal
#'
#' @importFrom glmmTMB fixef VarCorr sigma family_params
extractParams <- function(model, type, slope = FALSE) {
  fixed <- glmmTMB::fixef(model)$cond
  a <- fixed[1]
  b <- fixed[2]
  sigma2 <- attr(glmmTMB::VarCorr(model)$cond$strain, "stddev")^2
  phi <- glmmTMB::sigma(model)

  if (slope) {
    attrb <- attributes(glmmTMB::VarCorr(model)[[c("cond","strain")]])
    sig11 <- attrb$stddev[1]^2
    sig22 <- attrb$stddev[2]^2
    cor <- attrb$correlation[1,2]
    sig12 <- attrb$stddev[1] * attrb$stddev[2] * cor
  }

  if (type == "NB" && slope) {
    return(c(a, b, sig11, sig22, sig12, phi))
  } else if (type == "NB" && !slope) {
    return(c(a, b, sigma2, phi))
  } else if (type == "CP" && slope) {
    p <- glmmTMB::family_params(model)
    return(c(a, b, sig11, sig22, sig12, p, phi))
  } else if(type == "CP" && !slope) {
    p <- glmmTMB::family_params(model)
    return(c(a, b, sigma2, p, phi))
  }
}

#' Fit Generalized Linear Mixed-Effects Model to Count Matrix
#'
#' This function fits a generalized linear mixed-effects model to each row of a provided count matrix using the `glmmTMB` package.
#'
#' @param countMatrix A matrix where rows represent features and columns represent samples.
#' @param X A numeric vector representing the covariate.
#' @param type A character string specifying the type of model to fit. Either "NB" (negative binomial) or "CP" (tweedie).
#' @param slope Logical. If TRUE, the random effect will also have a slope component with the covariate. Default is FALSE.
#' @param parallel integer Set number of OpenMP threads to evaluate the negative log-likelihood in parallel
#'
#' @return A matrix with rows corresponding to features in `countMatrix` and columns representing the model parameters. Row names of the returned matrix match the row names of the input `countMatrix`.
#'
#' @seealso \code{\link[glmmTMB]{glmmTMB}}
#'
#' @keywords internal
fitGeneralizedModel <- function(countMatrix, X, type, slope=FALSE, parallel=1) {
  print(parallel)
  IDs <- rownames(countMatrix)
  rdfx <- as.factor(sapply(strsplit(colnames(countMatrix), "_"), `[`, 1))

  familyType <- if(type == "NB")  glmmTMB::nbinom2 else glmmTMB::tweedie

  fitModelForRow <- function(idx) {
    countVector <- countMatrix[idx, ]
    countData <- data.frame(expr = as.numeric(countVector),
                            covariate = X, strain = rdfx)
    formulaStr <- if(slope)
      stats::formula(expr ~ 1 + covariate + (1 + covariate | strain))
    else
      stats::formula(expr ~ 1 + covariate + (1 | strain))

    fullModel <- try({glmmTMB::glmmTMB(formula = formulaStr,
                                       data = countData,
                                       family = familyType,
                                       control = glmmTMB::glmmTMBControl(
                                         parallel=parallel,
                                         #optCtrl = list(iter.max=1e3,eval.max=1e3)
                                       ),
                                       REML = TRUE)}, silent = TRUE)
    if (!inherits(fullModel, "try-error")) {
      return(extractParams(fullModel, type, slope))
    } else {
      colLength <- ifelse(type == "NB", ifelse(slope, 6, 4), ifelse(slope, 7, 5))
      message(paste("Fitting problem for feature", idx, "returning NA"))
      return(rep(NA, colLength))
    }
  }

  # Apply the fitModelForRow function to each row
  params <- t(pbapply::pbsapply(1:nrow(countMatrix), fitModelForRow))

  # Define column names
  baseNames <- c("alpha", "beta")
  nbNames <- if(slope) c("sigma[11]^2", "sigma[22]^2", "sigma[12]", "phi")
  else c("sigma^2", "phi")
  cpNames <- if(slope) c("sigma[11]^2", "sigma[22]^2", "sigma[12]", "p", "phi")
  else c("sigma^2", "p", "phi")

  colNames <- c(baseNames, if(type == "NB") nbNames else cpNames)
  colnames(params) <- colNames
  rownames(params) <- IDs

  return(params)
}

#' Fit Negative Binomial Generalized Linear Mixed-Effects Model to Count Matrix
#'
#' Fits a negative binomial generalized linear mixed-effects model to each row of a provided count matrix using the `fitGeneralizedModel` function.
#'
#' @param countMatrix A matrix where rows represent features and columns represent samples.
#' @param X A numeric vector representing the covariate.
#' @param slope Logical. If TRUE, the random effect will also have a slope component with the covariate. Default is FALSE.
#' @param parallel Set number of OpenMP threads to evaluate the negative log-likelihood in parallel
#'
#' @return A matrix with rows corresponding to features in `countMatrix` and columns representing the model parameters. Row names of the returned matrix match the row names of the input `countMatrix`.
#' @export
#'
#' @seealso \code{\link{fitGeneralizedModel}}
fitNBmodel <- function(countMatrix, X, slope=FALSE, parallel=1)
  fitGeneralizedModel(countMatrix, X, "NB", slope=slope, parallel=parallel)

#' Fit Compound Poisson Generalized Linear Mixed-Effects Model to Count Matrix
#'
#' Fits a compound Poisson generalized linear mixed-effects model to each row of a provided count matrix using the `fitGeneralizedModel` function.
#'
#' @param countMatrix A matrix where rows represent features and columns represent samples.
#' @param X A numeric vector representing the covariate.
#' @param slope Logical. If TRUE, the random effect will also have a slope component with the covariate. Default is FALSE.
#' @param parallel Set number of OpenMP threads to evaluate the negative log-likelihood in parallel
#'
#' @return A matrix with rows corresponding to features in `countMatrix` and columns representing the model parameters. Row names of the returned matrix match the row names of the input `countMatrix`.
#' @export
#'
#' @seealso \code{\link{fitGeneralizedModel}}
#'
fitCPmodel <- function(countMatrix, X, slope=FALSE, parallel=1)
  fitGeneralizedModel(countMatrix, X, "CP", slope=slope, parallel=parallel)
