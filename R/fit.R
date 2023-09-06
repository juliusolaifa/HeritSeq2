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

fitGeneralizedModel <- function(countMatrix, X, type, slope=FALSE) {
  IDs <- rownames(countMatrix)
  rdfx <- as.factor(sapply(strsplit(colnames(countMatrix), "_"), `[`, 1))

  familyType <- if(type == "NB")  "nbinom2" else "tweedie"

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

fitNBmodel <- function(countMatrix, X, slope=FALSE)
  fitGeneralizedModel(countMatrix, X, "NB", slope=slope)

fitCPmodel <- function(countMatrix, X, slope=FALSE)
  fitGeneralizedModel(countMatrix, X, "CP", slope=slope)
