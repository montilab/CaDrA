
test_that("revealer_genescore_mat generates a list with one column", {
  
  mat <- matrix(c(1,0,1,0,0,0,0,0,1,0, 
                  0,0,1,0,1,0,1,0,0,0,
                  0,0,0,0,1,0,1,0,1,0), nrow=3)
  
  rownames(mat) <- c("TP_1", "TP_2", "TP_3")
  colnames(mat) <- 1 : ncol(mat)
  
  input_score = rnorm(n = ncol(mat))
  names(input_score) <- colnames(mat)
  
  
  result <- suppressWarnings(
    revealer_rowscore(
      FS = mat, 
      input_score = input_score, 
      seed_names = NULL, 
      assoc_metric = "IC"
    )
  )
  
  testthat::expect_equal(dim(result), c(3,1))
  testthat::expect_type(result, "double")

})

