test_that("Validation tests for generateCovariate function", {
  # Test for empty or NULL 'strainSizes'
  expect_error(generateCovariate(NULL), "The 'strainSizes' vector cannot be empty.")
  expect_error(generateCovariate(numeric(0)), "The 'strainSizes' vector cannot be empty.")

  # Test to ensure 'strainSizes' is a numeric vector
  expect_error(generateCovariate(c("a", "b")), "'strainSizes' must be a numeric vector.")

  # Test for negative values in 'strainSizes'
  expect_error(generateCovariate(c(-5, 7)))

  # Test for non-integer values in 'strainSizes'
  expect_error(generateCovariate(c(5.5, 7.2)))

  # Test for invalid distribution function
  expect_error(generateCovariate(c(5, 7), distribution="unknown_func"))

  # Test for invalid parameters for distribution
  expect_error(generateCovariate(c(5, 7), distribution="rnorm", monkey="banana"))
})

test_that("Functionality tests for generateCovariate function", {
  # Test default distribution
  sample_sizes <- c(10, 15, 20)
  covariate <- generateCovariate(sample_sizes)
  expect_true(all(covariate %in% c(0, 1)))
  expect_equal(length(covariate), sum(sample_sizes))

  # Test normal distribution
  sample_sizes <- c(5, 7)
  covariate <- generateCovariate(sample_sizes, distribution="rnorm")
  expect_true(is.numeric(covariate))
  expect_equal(length(covariate), sum(sample_sizes))

  # Test uniform distribution
  sample_sizes <- c(4, 6)
  min_val <- 0
  max_val <- 1
  covariate <- generateCovariate(sample_sizes, distribution="runif")
  expect_true(all(covariate >= min_val))
  expect_true(all(covariate <= max_val))
  expect_equal(length(covariate), sum(sample_sizes))

  # Test large sample sizes
  sample_sizes <- rep(1e6, 5) # 5 million total samples
  covariate <- generateCovariate(sample_sizes)
  expect_equal(length(covariate), sum(sample_sizes))
})

test_that("getPartOfColNames works correctly", {

  # Typical Cases
  testMatrix1 <- matrix(1:6, ncol=3)
  colnames(testMatrix1) <- c("name#1#test", "name#2#test", "name#3#test")
  expect_identical(getPartOfColNames(testMatrix1, 2, splitSymbol="#"), c("1", "2", "3"))

  # Edge Cases
  # Test when splitSymbol is not found
  testMatrix2 <- matrix(1:6, ncol=3)
  colnames(testMatrix2) <- c("name.1.test", "name.2.test", "name.3.test")
  expect_error(getPartOfColNames(testMatrix2, 2, splitSymbol="#"), "Error: splitSymbol not found in one or more column names.")

  # Test when column names are missing
  testMatrix3 <- matrix(1:6, ncol=3)
  expect_error(getPartOfColNames(testMatrix3, 2, splitSymbol="#"), "Error: countMatrix does not have column names.")

  # Test when partIndex is out of bounds
  testMatrix4 <- matrix(1:6, ncol=3)
  colnames(testMatrix4) <- c("name#1", "name#2", "name#3")
  expect_identical(is.na(getPartOfColNames(testMatrix4, 3, splitSymbol="#")), c(TRUE, TRUE, TRUE))

  # Potential Error Cases
  # Test with non-matrix or data.frame input
  expect_error(getPartOfColNames(c("a", "b", "c"), 2), "Error: countMatrix should be a matrix.")

  # Test with non-numeric partIndex
  expect_error(getPartOfColNames(testMatrix1, "two", splitSymbol="#"), "Error: partIndex should be a positive integer.")

  # Test with empty splitSymbol
  expect_error(getPartOfColNames(testMatrix1, 2, splitSymbol=""), "Error: splitSymbol should be a non-empty string.")

})
