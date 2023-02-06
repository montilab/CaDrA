
test_that("revealer_rowscore returns correct results", {
  
  
  mat <- matrix(c(1,0,1,0,0,0,0,0,1,0, 
                  0,0,1,0,1,0,1,0,0,0,
                  0,0,0,0,1,0,1,0,1,0), nrow=3)
  
  colnames(mat) <- 1:10
  row.names(mat) <- c("TP_1", "TP_2", "TP_3")
  
  set.seed(42)
  input_score = rnorm(n = ncol(mat))
  names(input_score) <- colnames(mat)
  
  result <- revealer_rowscore(
    FS_mat = mat, 
    input_score = input_score, 
    assoc_metric = "IC"
  )
  
  testthat::expect_identical(dim(result), c(3L,1L))
  testthat::expect_type(result, "double")
  testthat::expect_identical(rownames(result), row.names(mat))
  testthat::expect_identical(colnames(result), "score")
  testthat::expect_equal(round(result[,1], 3), c(TP_1=0.457, TP_2=-0.427, TP_3=0.387) )
  
  result <- revealer_rowscore(
    FS_mat = mat, 
    input_score = input_score, 
    assoc_metric = "COR"
  )
  testthat::expect_identical(dim(result), c(3L,1L))
  testthat::expect_type(result, "double")
  testthat::expect_identical(rownames(result), row.names(mat))
  testthat::expect_identical(colnames(result), "score")
  testthat::expect_equal(round(result[,1], 3), c(TP_1=0.592, TP_2=-0.399, TP_3=0.507) )
  
  
})

## --------------------------------------------------- ##
test_that("revealer_rowscore issues error messages when needed", {
  
  FS_mat <-  data.frame(a = rnorm(10), b = rnorm (10) )
  set.seed(42)
  input_score <- rnorm(n = ncol(FS_mat))
  names(input_score) <- colnames(FS_mat)
  
  expect_error( revealer_rowscore(
    FS_mat = FS_mat,  
    input_score = input_score)
  )
  
  
  mat <- matrix(c(1,0,1,0,0,0,0,0,1,0, 
                  0,0,1,0,1,0,1,0,0,0,
                  0,0,0,0,1,0,1,0,1,0), nrow=3)
  
  colnames(mat) <- paste0("COL_", 1:10)
  row.names(mat) <- c("TP_1", "TP_2", "TP_3")
  
  set.seed(42)
  input_score <- rnorm(n = ncol(mat))
  names(input_score) <- colnames(mat)
  seed_names <- as.character(1:10)
  expect_error( revealer_rowscore(
    FS_mat = mat, 
    input_score = input_score, 
    seed_names = seed_names,
    assoc_metric = "IC")
  )
  
  
  
  
})

## --------------------------------------------------- ##
test_that("revealer_score returns correct results", {
  
  set.seed(42)
  input_score <- rnorm(n = 10)
  x <- c(1,0,1,0,0,0,0,0,1,0)
  seed_vector <- rep(0, 10)
    
  result <- revealer_score(
    x = input_score,
    y = x,
    z = seed_vector,
    assoc_metric = "IC")
  
  testthat::expect_type(result, "double")
  testthat::expect_identical(names(result), "score")
  testthat::expect_equal(round(result, 3), c(score=0.44) )
  
  
})

## --------------------------------------------------- ##
test_that("cond_assoc returns correct results", {
  
  set.seed(42)

  result <- suppressWarnings(cond_assoc(
    x = rnorm(n = 10),
    y = c(1,0,1,0,0,0,0,0,1,0),
    z = rep(0, 10),
    metric = "IC"))
  
  testthat::expect_type(result, "double")
  testthat::expect_equal(round(result, 3), c(0.44) )
  
  
})
  

## --------------------------------------------------- ##
test_that("mutual_inf_v2 returns correct results", {
  
  set.seed(42)
  
  result <- suppressWarnings(mutual_inf_v2(
    x = rnorm(n = 10),
    y = x + rnorm(n = 10, 0, 0.01)))
  
  testthat::expect_type(result, "list")
  testthat::expect_identical(names(result), 
                             c("MI", "SMI", "HXY", "HX", "HY", "NMI", "IC") )
  
})



## --------------------------------------------------- ##
test_that("cond_assoc issues appropriate error messages", {
  
  
  
})
