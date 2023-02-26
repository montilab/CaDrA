
test_that("wilcox_rowscore returns correct results", {
  
  
  mat <- matrix(c(1,0,1,0,0,0,0,0,1,0, 
                  0,0,1,0,1,0,1,0,0,0,
                  0,0,0,0,1,0,1,0,1,0), nrow=3)
  
  colnames(mat) <- 1:10
  row.names(mat) <- c("TP_1", "TP_2", "TP_3")
  
  set.seed(42)
  input_score = rnorm(n = ncol(mat))
  names(input_score) <- colnames(mat)
  
  result <- wilcox_rowscore(
    FS_mat = mat, 
    input_score = input_score, 
    alternative = "less"
  )
  
  testthat::expect_length(result, 3L)
  testthat::expect_type(result, "double")
  testthat::expect_identical(result, c(TP_2=13,TP_3=5,TP_1=3))
  
})

## --------------------------------------------------- ##
test_that("wilcox_rowscore issues error messages when needed", {
  
  FS_mat <-  data.frame(a = rnorm(10), b = rnorm (10) )
  set.seed(42)
  input_score = rnorm(n = ncol(FS_mat))
  names(input_score) <- colnames(FS_mat)
  
  expect_error( wilcox_rowscore(
    FS_mat = FS_mat,  
    input_score = input_score)
  )
  
  
})


## --------------------------------------------------- ##
test_that("wilcox_score issues correct errors", {
  
  row <- c(1, 0, 1, 0, 0, 0, 0, 0, 1, 0)
  ranks <- 1:length(row)
  x = ranks[which(row == 1)]
  y = ranks[which(row == 0)]
  
  mu = NA
  expect_error( wilcox_score(x, y, mu) )
  
  mu = 0
  expect_error( wilcox_score(x=c(), y, mu) )
  expect_error( wilcox_score(x=x, y=c(), mu) )
  
  
})
  
## --------------------------------------------------- ##
test_that("wilcox_score returns correct results", {
  
  row <- c(1, 0, 1, 0, 0, 0, 0, 0, 1, 0)
  ranks <- 1:length(row)
  x = ranks[which(row == 1)]
  y = ranks[which(row == 0)]
  mu = 0
  
  result <- wilcox_score(x=x, y=y, mu=mu)
  expect_named(result, c("score.W", "p_value" ))
  expect_equal( result[1], c("score.W"=7.0) )
  
  expect_equal( result[2], 
                c("p_value"= pnorm( (7 - 21/2 + 0.5)/sqrt(21 * 11/12)) ) )
  
})
