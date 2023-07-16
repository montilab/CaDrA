
test_that("ks_rowscore returns correct results", {
  
  mat <- matrix(c(1,0,1,0,0,0,0,0,1,0, 
                  0,0,1,0,1,0,1,0,0,0,
                  0,0,0,0,1,0,1,0,1,0), nrow=3)
  
  colnames(mat) <- 1:10
  row.names(mat) <- c("TP_1", "TP_2", "TP_3")
  
  set.seed(42)
  input_score = rnorm(n = ncol(mat))
  names(input_score) <- colnames(mat)
  
  result <- ks_rowscore(
    FS = mat, 
    input_score = input_score, 
    meta_feature = NULL,
    alternative = "less"
  )
  
  testthat::expect_length(result, 3L)
  testthat::expect_type(result, "double")
  testthat::expect_equal(result, c(TP_1=0.5,TP_2=-0.6,TP_3=0.4))
  
  weights <- c(0.5, 0.75, 0.25, 1, 0.5, 0.25, 0.75, 0.5, 0.5, 0 )
  names(weights) <- colnames(mat)
  result <- ks_rowscore(
    FS = mat,  
    input_score = input_score,
    meta_feature = NULL,
    weights = weights
  )
  
  testthat::expect_length(result, 3L)
  testthat::expect_type(result, "double")
  testthat::expect_named(result, c("TP_1", "TP_2", "TP_3"))
  testthat::expect_equal(round(result,7), 
                         c(TP_1=0.7142857, TP_2=0.8750000, TP_3=0.6666667))
  
})

## --------------------------------------------------- ##
test_that("ks_rowscore issues error messages when needed", {
  
  FS <-  data.frame(a = rnorm(10), b = rnorm (10) )
  set.seed(42)
  input_score <- rnorm(n = ncol(FS))
  names(input_score) <- colnames(FS)
  
  expect_error( 
    ks_rowscore(
      FS = FS,  
      input_score = input_score,
      meta_feature = NULL,
    )
  )
  
  mat <- matrix(c(1,0,1,0,0,0,0,0,1,0, 
                  0,0,1,0,1,0,1,0,0,0,
                  0,0,0,0,1,0,1,0,1,0), nrow=3)
  
  colnames(mat) <- paste0("COL_", 1:10)
  
  set.seed(42)
  input_score <- rnorm(n = ncol(mat))
  names(input_score) <- colnames(mat)
  weights <- runif(n = ncol(mat))
  
  expect_error( 
    ks_rowscore(
      FS = mat,  
      input_score = input_score,
      meta_feature = NULL,
      weights = weights
    )
  )
  
})

