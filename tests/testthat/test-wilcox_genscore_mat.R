test_that(
  "wilcox_genescore_mat generates results consistent with wilcox() function ", {
    
    x <- c(15, 20, 39, 42, 44)
    y <- (1:100)[-x]
    res.wilcox<- wilcox_score(x=x, y=y, alternative="less")
    
    CORRECTION=-0.5
    
    n.x <- 5
    n.y <- 95
    r <- c(x,y)
    W <- sum(r[seq_along(x)]) - n.x * (n.x + 1)/2
    SIGMA <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) ))
    z <- W - n.x * n.y/2
    
    z <- (z - CORRECTION)/SIGMA
    p <- pnorm(z)
    
    
    testthat::expect_equal(unname(res.wilcox["score.W"]), W)
    testthat::expect_equal(unname(res.wilcox["p_value"]), p)
    
    
  })

test_that("wilcox_genescore_mat generates a matrix with 2 rows", {
  
  mat <- matrix(c(1,0,1,0,0,0,0,0,1,0, 
                  0,0,1,0,1,0,1,0,0,0,
                  0,0,0,0,1,0,1,0,1,0), nrow=3)
  
  colnames(mat) <- 1:10
  row.names(mat) <- c("TP_1", "TP_2", "TP_3")
  
  input_score = rnorm(n = ncol(mat))
  names(input_score) <- colnames(mat)
  
  result<- wilcox_rowscore(mat, input_score, alt="less")
  
  testthat::expect_equal(dim(result), c(3,2))
  testthat::expect_type(result, "double")

  
})

