
test_that("revealer_rowscore returns correct results", {
  
  mat <- matrix(c(1,0,1,0,0,0,0,0,1,0, 
                  0,0,1,0,1,0,1,0,0,0,
                  0,0,0,0,1,0,1,0,1,0), nrow=3)
  
  colnames(mat) <- 1:10
  row.names(mat) <- c("TP_1", "TP_2", "TP_3")
  
  # Set seed
  set.seed(42)
  
  input_score = rnorm(n = ncol(mat))
  names(input_score) <- colnames(mat)
  
  result <- revealer_rowscore(
    FS = mat, 
    input_score = input_score, 
    meta_feature = NULL,
    assoc_metric = "IC"
  )
  
  testthat::expect_length(result, 3L)
  testthat::expect_type(result, "double")
  testthat::expect_identical(sort(names(result)), sort(row.names(mat)))
  testthat::expect_equal(round(result, 7), 
                         c(TP_1=0.4573207, TP_2=-0.4266714, TP_3=0.3867053))
  
  result <- revealer_rowscore(
    FS = mat, 
    input_score = input_score, 
    meta_feature = NULL,
    assoc_metric = "COR"
  )
  
  testthat::expect_length(result, 3L)
  testthat::expect_type(result, "double")
  testthat::expect_identical(sort(names(result)), sort(row.names(mat)))
  testthat::expect_equal(round(result, 7), 
                         c(TP_1=0.5924383, TP_2=-0.3985206, TP_3=0.5067436))
  
})

## --------------------------------------------------- ##
test_that("revealer_score returns correct results", {
  
  # Set seed
  set.seed(42)
  
  input_score <- rnorm(n = 10)
  x <- c(1,0,1,0,0,0,0,0,1,0)
  meta_vector <- rep(0, 10)
    
  result <- revealer_score(
    x = input_score,
    y = x,
    z = meta_vector,
    assoc_metric = "IC"
  )
  
  testthat::expect_type(result, "double")
  testthat::expect_named(result, "score")
  testthat::expect_equal(round(result, 7), c(score=0.4400776) )
  
})

## --------------------------------------------------- ##
test_that("cond_assoc returns correct results", {
  
  # Set seed
  set.seed(42)

  result <- suppressWarnings(
    cond_assoc(
      x = rnorm(n = 10),
      y = c(1,0,1,0,0,0,0,0,1,0),
      z = rep(0, 10),
      metric = "IC"
    )
  )
  
  testthat::expect_type(result, "double")
  testthat::expect_equal(round(result, 7), c(0.4400776) )
  
})
  

## --------------------------------------------------- ##
test_that("mutual_inf_v2 returns correct results", {
  
  # Set seed
  set.seed(42)
  
  x <- rnorm(n = 10)
  y <- x + rnorm(n = 10, 0, 0.01)
  
  result <- suppressWarnings(
    mutual_inf_v2(
      x = x,
      y = y
    )
  )
  
  testthat::expect_type(result, "list")
  testthat::expect_named(result, 
                         c("MI", "SMI", "HXY", "HX", "HY", "NMI", "IC") )
  
  testthat::expect_equal(round(result$MI, 5), 1.11879 )
  testthat::expect_equal(round(result$SMI, 5), 1.11879 )
  testthat::expect_equal(round(result$HXY, 5), 0.42585 )
  testthat::expect_equal(round(result$HX, 6), 0.779993 )
  testthat::expect_equal(round(result$HY, 7), 0.7646444 )
  testthat::expect_equal(round(result$NMI, 6), 2.627207 )
  testthat::expect_equal(round(result$IC, 6),  0.945137 )
  
})


