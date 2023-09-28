
test_that("custom_rowscore returns correct results", {
  
  mat <- matrix(c(1,0,1,0,0,0,0,0,1,0, 
                  0,0,1,0,1,0,1,0,0,0,
                  0,0,0,0,1,0,1,0,1,0), nrow=3)
  
  colnames(mat) <- 1:10
  row.names(mat) <- c("TP_1", "TP_2", "TP_3")
  
  set.seed(42)
  input_score = rnorm(n = ncol(mat))
  names(input_score) <- colnames(mat)
  custom_function <- wilcox_rowscore
  
  result <- custom_rowscore(
    FS = mat, 
    input_score = input_score,
    meta_feature = NULL,
    custom_function = wilcox_rowscore,
    custom_parameters = list(alternative="less")
  )
  
  testthat::expect_length(result, 3L)
  testthat::expect_type(result, "double")
  testthat::expect_identical(result, c(TP_1=3, TP_2=13, TP_3=5))
  
})

## --------------------------------------------------- ##
test_that("custom_rowscore issues error messages when needed", {
  
  mat <- matrix(c(1,0,1,0,0,0,0,0,1,0, 
                  0,0,1,0,1,0,1,0,0,0,
                  0,0,0,0,1,0,1,0,1,0), nrow=3)
  
  colnames(mat) <- 1:10
  row.names(mat) <- c("TP_1", "TP_2", "TP_3")
  
  set.seed(42)
  input_score = rnorm(n = ncol(mat))
  names(input_score) <- colnames(mat)
  custom_function <- wilcox_rowscore

  # error if custom function is not given 
  expect_error( 
    custom_rowscore(
      FS = mat, 
      input_score = input_score,
      meta_feature = NULL,
      custom_function = "wilcox_rowscore",
      method = "custom_pval", 
      custom_parameters = list(alternative="less") 
    )
  )
  
  # error if custom parameters are not given in a list
  expect_error( 
    custom_rowscore(
      FS = mat, 
      input_score = input_score,
      meta_feature = NULL,
      custom_function = custom_function,
      method = "custom_pval", 
      custom_parameters = c(alternative="less")
    )
  )
  
  # error if function does not have "input_score" and "meta_feature" arguments
  err_function <- function(FS){ return(0) }
  expect_error( 
    custom_rowscore(
      FS = mat, 
      input_score = input_score,
      meta_feature = NULL,
      custom_function = err_function,
      method = "custom_pval"
    )
  )
  
  # error if function does not have "FS" and "meta_feature" arguments
  err_function <- function(input_score){ return(0) }
  expect_error( 
    custom_rowscore(
      FS = mat, 
      input_score = input_score,
      meta_feature = NULL,
      custom_function = err_function,
      method = "custom_pval"
    )
  )
  
})
