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

#' Fit Generalized Model to Each Row of Count Matrix
#'
#' This function fits a generalized model (either Negative Binomial or Tweedie)
#' to each row of a given count matrix. It returns the model results and
#' optionally captures any warnings that are generated during the fitting process.
#'
#' @param countMatrix A matrix where rows represent features and columns represent samples.
#' @param X A numeric vector representing the covariate.
#' @param slope Logical. Determines whether the random effect should include a slope component with the covariate. Default is FALSE.
#' @param type Character string indicating the model type: "NB" for Negative Binomial or "CP" for Compound Poisson (Tweedie).
#' @param parallel Numeric value indicating the number of OpenMP threads to use for parallel processing.
#' @param returnWarnings Logical. If TRUE, the function will return any warnings that were generated during the model fitting process.
#'
#' @return A list where each element corresponds to the result of fitting the model
#' to a row of the countMatrix. If \code{returnWarnings} is TRUE, this list will contain
#' any warnings that were generated. Otherwise, it will contain model parameters or NA values for failed fits.
#'
#' @examples
#' # This is a placeholder. Update with a relevant example:
#' # countMatrix <- your_example_data
#' # X <- your_example_data
#' # fitModelCommon(countMatrix, X, slope=TRUE, type="NB", parallel=2, returnWarnings=TRUE)
#'
#' @seealso
#' \code{\link{captureModelWarnings}}, \code{\link{fitGeneralizedModel}}
#' @export
fitModelCommon <- function(countMatrix, X, slope, type, parallel, returnWarnings = FALSE) {
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

    # internal function to fit the model
    fit_model <- function(formulaStr, countData, familyType, parallel) {
      glmmTMB::glmmTMB(formula = formulaStr,
                       data = countData,
                       family = familyType,
                       control = glmmTMB::glmmTMBControl(optCtrl =
                                                           list(iter.max=300,
                                                                eval.max=400),
                                                         parallel=parallel),
                       REML = TRUE)}

    # handle warnings or model
    if (returnWarnings) {
      model <- tryCatch({
         model_object <- fit_model(formulaStr, countData, familyType, parallel)
        return(paste(idx, "No Warning"))
      }, warning = function(w) {
        return(paste(idx, w$message))
      },
      silent = TRUE)
    }
    else {
      fit_time <- system.time(
        model <- try({
          fit_model(formulaStr, countData, familyType, parallel)
        },silent = TRUE)
      )
      time_taken <- fit_time["elapsed"]

      # Print the time taken for this row
      cat(paste("Gene", idx, "- Time taken for glmmTMB fit:", round(time_taken, 4), "seconds\n"))

      if (!inherits(model, "try-error")) {
        return(extractParams(model, type, slope))
      } else {
        colLength <- ifelse(type == "NB", ifelse(slope, 6, 4), ifelse(slope, 7, 5))
        return(rep(NA, colLength))
      }
    }
  }

  #results <- microbenchmark::microbenchmark()
  return(pbapply::pblapply(1:nrow(countMatrix), fitModelForRow))
}


#' Capture Model Warnings
#'
#' This function captures warnings generated during the model fitting process
#' and returns them in a structured format.
#'
#' @param countMatrix A matrix where rows represent features and columns represent samples.
#' @param X A numeric vector representing the covariate.
#' @param type Character indicating the model type (e.g., "NB" for Negative Binomial).
#' @param slope Logical indicating whether the random effect will have a slope component with the covariate. Default is FALSE.
#' @param parallel Numeric. Specifies the number of OpenMP threads to use for parallel processing.
#'
#' @return A dataframe with the captured warnings.
#' Each row corresponds to a feature in the countMatrix and contains the corresponding warning message.
#'
#' @examples
#' # This is a placeholder. Update with a relevant example:
#' # countMatrix <- your_example_data
#' # X <- your_example_data
#' # captureModelWarnings(countMatrix, X, type="NB")
#'
#' @seealso
#' \code{\link{fitModelCommon}}, \code{\link{fitGeneralizedModel}}
#' @export
captureModelWarnings <- function(countMatrix, X, type, slope=FALSE, parallel=1) {
  warnings_list <- fitModelCommon(countMatrix, X, slope, type, parallel, returnWarnings=TRUE)
  # Extract only the unique warnings for each list element
  warnings_list <- lapply(warnings_list, function(warn) unique(warn))
  df <- data.frame(column=unlist(warnings_list))

  # Extract the initial integer and the rest of the string
  df$integer <- stringr::str_extract(df$column, "^[0-9]+")
  df$warning <- sub("^[0-9]+ ", "", df$column)

  rownames(df) <- df$integer
  df$column <- NULL
  df$integer <- NULL

  return(df)
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

#' @return A matrix with rows corresponding to features in `countMatrix` and columns representing the model parameters. Row names of the returned matrix match the row names of the input `countMatrix`.
#'
#' @seealso \code{\link[glmmTMB]{glmmTMB}}
#'
#' @keywords internal
fitGeneralizedModel <- function(countMatrix, X, type, slope = FALSE, parallel = 1) {
  paramsList <- fitModelCommon(countMatrix, X, slope, type, parallel, returnWarnings = FALSE)

  # Transforming the list into a matrix
  params <- do.call(rbind, paramsList)

  # Define column names
  baseNames <- c("alpha", "beta")
  nbNames <- if(slope) c("sigma[11]^2", "sigma[22]^2", "sigma[12]", "phi")
  else c("sigma^2", "phi")
  cpNames <- if(slope) c("sigma[11]^2", "sigma[22]^2", "sigma[12]", "p", "phi")
  else c("sigma^2", "p", "phi")

  colNames <- c(baseNames, if(type == "NB") nbNames else cpNames)
  colnames(params) <- colNames
  rownames(params) <- rownames(countMatrix)

  return(params)
}

# This function handles whether  fitNBModel and fitCPModel should only fit model
# or also capture warnings
.fitModelHelper <- function(countMatrix, X, modelType, slope = FALSE, parallel = 1, reportWarning = FALSE) {

  if(modelType == "CP" && parallel == 1) {
    message("Consider training Compound Poisson model in parallel -- set : parallel > 1")
  }

  suppressWarnings({
    fit <- fitGeneralizedModel(countMatrix, X, modelType, slope = slope, parallel = parallel)
    output <- list(params = fit)
    if(reportWarning) {
      warn <- captureModelWarnings(countMatrix, X, modelType, slope = FALSE, parallel = 1)
      output$warnings <- warn
    }
  })

  return(output)
}


#' Fit a Negative Binomial Model to a Count Matrix
#'
#' Fits a Negative Binomial model to each row of a provided count matrix.
#'
#' @param countMatrix A matrix where rows represent features and columns represent samples.
#' @param X A numeric vector representing the covariate.
#' @param slope Logical. If TRUE, the random effect will also have a slope component
#'        with the covariate. Default is FALSE.
#' @param parallel Set number of OpenMP threads to evaluate the negative log-likelihood in parallel.
#' @param reportWarning Logical. If TRUE then the model is fit twice to obtain the warning logs.
#'
#' @return A list containing the model parameters and, if \code{reportWarning} is TRUE,
#'         the warnings generated during model fitting.
#'
#' @seealso \code{\link{fitGeneralizedModel}}, \code{\link{captureModelWarnings}}
#' @export
fitNBmodel <- function(countMatrix, X, slope = FALSE, parallel = 1, reportWarning = FALSE) {
  print("fitting with NB")
  .fitModelHelper(countMatrix, X, "NB", slope, parallel, reportWarning)
}


#' Fit a Compound Poisson Model to a Count Matrix
#'
#' Fits a Compound Poisson model to each row of a provided count matrix.
#'
#' @param countMatrix A matrix where rows represent features and columns represent samples.
#' @param X A numeric vector representing the covariate.
#' @param slope Logical. If TRUE, the random effect will also have a slope component
#'        with the covariate. Default is FALSE.
#' @param parallel Set number of OpenMP threads to evaluate the negative log-likelihood in parallel.
#' @param reportWarning Logical. If TRUE then the model is fit twice to obtain the warning logs.
#'
#' @return A list containing the model parameters and, if \code{reportWarning} is TRUE,
#'         the warnings generated during model fitting.
#'
#' @seealso \code{\link{fitGeneralizedModel}}, \code{\link{captureModelWarnings}}
#' @export
fitCPmodel <- function(countMatrix, X, slope = FALSE, parallel = 1, reportWarning = FALSE) {
  print("fitting with CP")
  .fitModelHelper(countMatrix, X, "CP", slope, parallel, reportWarning)
}
