#' Generate binary covariate

#' @param strainSizes A numeric vector of integers.
#'  Each element denotes the number of samples within a corresponding strain.
#'  approximately equal numbers of 0's and 1's are generated within each strain.

#'  @noRd
rgen01 <- function(strainSizes) {
  unlist(lapply(strainSizes, function(x) sample(rep_len(c(0,1),length.out=x))))
}

#' Covariate Generation from Specified Distribution

#' @param strainSizes A numeric vector of integers.
#'  Each element denotes the number of samples within a corresponding strain.
#'  The length of this vector indicates the total number of strains.

#' @param distribution A function representing the desired distribution from
#'  which the covariate values will be generated. By default, the function
#'  generates a binary covariate. For functions in `R` use the generating functions
#'  e.g `rnorm`, `rpois`, `rgamma` e.t.c

#' @param ... Optional named parameters specific to the provided `distribution`
#' function

#' @return A numeric vector representing the generated covariate from the
#' specified distribution.
#' @export

#' @examples
#' strainSizes <- c(10, 15, 20)
#' generateCovariate(strainSizes)
#' generateCovariate(strainSizes, distribution = "rnorm")
#' generateCovariate(strainSizes, distribution = "rgamma", shape=2)
generateCovariate <- function(strainSizes, distribution="rgen01", ...) {
  # Check if 'strainSizes' is NULL or empty; if so, stop execution
  if(is.null(strainSizes) || length(strainSizes) == 0) {
    stop("The 'strainSizes' vector cannot be empty.")
  }

  # Check if 'strainSizes' is a numeric vector; if not, stop execution
  if(!is.numeric(strainSizes)) {
    stop("'strainSizes' must be a numeric vector.")
  }

  # Check for negative values in 'strainSizes'; if found, stop execution
  if(any(strainSizes < 0)) {
    stop("Elements in 'strainSizes' cannot be negative.")
  }

  # Ensure all elements in 'strainSizes' are integers; if not, stop execution
  if(any(floor(strainSizes) != strainSizes)) {
    stop("Elements in 'strainSizes' must be integers.")
  }

  # If the distribution is not 'rgen01', sum up the 'strainSizes'
  if (distribution != "rgen01") strainSizes = sum(strainSizes)

  # Try to generate values using the specified distribution
  # If there's an error, catch it and provide a descriptive message
  result <- tryCatch({
    do.call(distribution, list(strainSizes, ...))
  }, error = function(e) {
    stop(paste("Error in the distribution function:", e$message))
  })

  # Return the generated values
  return(result)
}


#' Extract Specific Parts of Column Names Based on a Delimiter
#'
#' This function takes a matrix, splits its column names based on a specified delimiter,
#' and then returns the part of the split name corresponding to a specified position.
#'
#' @param countMatrix A matrix containing the data. This matrix should have column names.
#' @param partIndex An integer specifying which part of the split column names to return.
#'        For example, if `partIndex` is 2 and the column name is "A_B_C", the function will return "B".
#' @param splitSymbol A string specifying the delimiter to use when splitting column names.
#'        Default is "_".
#'
#' @return A character vector containing the extracted parts of the column names.
#' @export
#'
#' @examples
#' testMatrix <- matrix(1:6, ncol=3)
#' colnames(testMatrix) <- c("name#1#test", "name#2", "name#3#extra#info")
#' getPartOfColNames(testMatrix, 2, splitSymbol="#")
getPartOfColNames <- function(countMatrix, partIndex, splitSymbol="_") {
  # Check for valid input types
  if (!is.matrix(countMatrix)) {
    stop("Error: countMatrix should be a matrix.")
  }
  if (!is.numeric(partIndex) || partIndex <= 0 || partIndex != round(partIndex)) {
    stop("Error: partIndex should be a positive integer.")
  }
  if (!is.character(splitSymbol) || nchar(splitSymbol) == 0) {
    stop("Error: splitSymbol should be a non-empty string.")
  }

  colNames <- colnames(countMatrix)
  if (is.null(colNames)) {
    stop("Error: countMatrix does not have column names.")
  }

  # splitting the string as many times as possible.
    splitNames <- stringr::str_split(colNames, splitSymbol)

    # Check if splitSymbol is not found in any of the column names
    if (any(sapply(splitNames, length) == 1)) {
      stop("Error: splitSymbol not found in one or more column names.")
    }

  # Extract the desired column (partIndex) from the matrix
    extractedParts <- sapply(splitNames, function(x) {
      if (length(x) >= partIndex) {
        return(x[partIndex])
      } else {
        return(NA)
      }
    })

    return(extractedParts)
}

#' Extract a Specific Part of Column Names and Convert to Factor
#'
#' This function takes a matrix, splits its column names based on a specified delimiter,
#' and then returns a factor representing the part of the split name corresponding to a specified position.
#'
#' @param countMatrix A matrix or data frame with column names.
#' @param partIndex An integer specifying which part of the split column names to return.
#' @param splitSymbol A string specifying the delimiter to use when splitting column names. Default is "_".
#'
#' @return A factor representing the extracted part of the column names.
#' @export
#' @examples
#' testMatrix <- matrix(1:6, ncol=3)
#' colnames(testMatrix) <- c("strainA_part1_result", "strainB_part2_result", "strainC_part3_result")
#' extractStrainFactor(testMatrix, 1)
extractStrainFactor <- function(countMatrix, partIndex, splitSymbol="_") {
  as.factor(getPartOfColNames(countMatrix, partIndex, splitSymbol))
}

#' Convert Specific Group of Column Names to Binary
#'
#' This function extracts a specified group from column names based on a delimiter, then
#' checks if each extracted group matches the target symbol(s). The result is a binary
#' representation, where matches are represented as `1` (or `trueValue`), and non-matches
#' as `0` (or `falseValue`).
#'
#' @param countMatrix A matrix or data frame with column names.
#' @param partIndex A positive integer specifying the position of the part in the split column names.
#' @param splitSymbol A character string specifying the delimiter to split the column names. Default is "_".
#' @param targetSymbols A character or character vector specifying the target group(s) to match. Default is "E".
#' @param trueValue Value to be returned for matching groups. Default is `1`.
#' @param falseValue Value to be returned for non-matching groups. Default is `0`.
#'
#' @return A numeric vector of binary values representing whether each column name's
#' specified part matches the target symbol(s).
#' @export
#' @examples
#' testMatrix <- matrix(1:6, ncol=3)
#' colnames(testMatrix) <- c("name_E_1", "name_F_2", "name_E_3")
#' extractGroupAsBinary(testMatrix, 2) # Expected: 1 0 1
extractGroupAsBinary <- function(countMatrix, partIndex, splitSymbol="_", targetSymbols ="E", trueValue=1, falseValue=0) {
  extractedGroups <- getPartOfColNames(countMatrix, partIndex, splitSymbol)
  ifelse(extractedGroups %in% targetSymbols, trueValue, falseValue)
}

#' Generate Histogram Plots for Parameters in a Matrix
#'
#' @param paramMatrix A matrix with numeric data.
#' @param col Color of the histogram bars. Default is "lightblue".
#' @param border Border color of the histogram bars. Default is "black".
#' @export
#'
#' @examples
#' set.seed(123) # Setting a seed for reproducibility
#' sampleMatrix <- matrix(rnorm(300), ncol=3)
#' colnames(sampleMatrix) <- c("Parameter_1", "Parameter_2", "Parameter_3")
#'
#' # Use the generateParamPlots function to create histograms
#' generateParamPlots(sampleMatrix)
generateParamPlots <- function(paramMatrix, col="lightblue", border="black") {
  # Error checks
  if (!is.matrix(paramMatrix)) {
    stop("paramMatrix should be a matrix.")
  }

  nParams <- ncol(paramMatrix)
  columnNames <- colnames(paramMatrix)

  # Calculate the grid for plotting
  nRows <- ceiling(nParams/2)
  op <- par(mfrow=c(nRows, 2), mar=c(4,4,2,1)+0.1) # Store original plotting parameters

  for (i in 1:nParams) {
    hist(paramMatrix[, i], main = parse(text = columnNames[i]),
         xlab = "Values", col = col, border = border)
  }

  # Resetting the plotting parameters to original values
  par(op)
}
