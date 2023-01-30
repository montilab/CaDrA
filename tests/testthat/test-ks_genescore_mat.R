test_that(
  "ks_genescore_mat generates results consistent with ks.test() function", {
    
    res.ks.test <- suppressWarnings(
      ks.test(
        x=1:10, c(1,3,9), 
        alternative="less", 
        exact=FALSE
      )
    )
    
    mat <- matrix(c(1,0,1,0,0,0,0,0,1,0), nrow=1)
    colnames(mat) <- 1:ncol(mat)
    rownames(mat) <- "f1"
    
    input_score <- ncol(mat):1
    names(input_score) <- 1:ncol(mat)
    
    res.cadra <- ks_rowscore(
      FS_mat=mat, 
      input_score=input_score, 
      weight=NULL, 
      alt="less"
    )
    
    testthat::expect_equal(res.cadra[1,1], unname(res.ks.test$statistic))
    testthat::expect_equal(res.cadra[1,2], unname(res.ks.test$p.value))
    
    
  })

test_that("ks_genescore_mat generates a matrix with 2 rows", {
  
  mat <- matrix(c(1,0,1,0,0,0,0,0,1,0, 
                  0,0,1,0,1,0,1,0,0,0,
                  0,0,0,0,1,0,1,0,1,0), nrow=3)
  
  colnames(mat) <- 1:ncol(mat)
  rownames(mat) <- c("f1", "f2", "f3")
  
  input_score <- ncol(mat):1
  names(input_score) <- 1:ncol(mat)
  
  result<- ks_rowscore(
    FS_mat=mat, 
    input_score=input_score, 
    weight=NULL, 
    alt="less"
  )
  
  testthat::expect_equal(dim(result), c(3,2))
  
  
})

