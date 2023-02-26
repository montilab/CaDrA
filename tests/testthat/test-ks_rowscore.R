
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
    FS_mat = mat, 
    input_score = input_score, 
    alternative = "less"
  )
  
  testthat::expect_length(result, 3L)
  testthat::expect_type(result, "double")
  testthat::expect_equal(result, c(TP_1=0.5,TP_3=0.4,TP_2=-0.6))
  
  weight <- c(0.5, 0.75, 0.25, 1, 0.5, 0.25, 0.75, 0.5, 0.5, 0 )
  names(weight) <- colnames(mat)
  result <- ks_rowscore(
    FS_mat = mat,  
    input_score = input_score,
    weight = weight)
  
  testthat::expect_length(result, 3L)
  testthat::expect_type(result, "double")
  testthat::expect_named(result, c("TP_2", "TP_1", "TP_3"))
  testthat::expect_equal(round(result,7), 
                         c(TP_2=0.8750000, TP_1=0.7142857, TP_3=0.6666667))
  
  
})

## --------------------------------------------------- ##
test_that("ks_rowscore issues error messages when needed", {
  
  FS_mat <-  data.frame(a = rnorm(10), b = rnorm (10) )
  set.seed(42)
  input_score <- rnorm(n = ncol(FS_mat))
  names(input_score) <- colnames(FS_mat)
  
  expect_error( ks_rowscore(
    FS_mat = FS_mat,  
    input_score = input_score)
  )
  
  
  mat <- matrix(c(1,0,1,0,0,0,0,0,1,0, 
                  0,0,1,0,1,0,1,0,0,0,
                  0,0,0,0,1,0,1,0,1,0), nrow=3)
  
  colnames(mat) <- paste0("COL_", 1:10)
  
  set.seed(42)
  input_score <- rnorm(n = ncol(mat))
  names(input_score) <- colnames(mat)
  weight <- runif(n = ncol(mat))
  
  expect_error( ks_rowscore(
    FS_mat = mat,  
    input_score = input_score,
    weight = weight)
  )
  
  
})

