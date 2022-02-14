test_that("ks_genescore_mat receives correct input arguments ", {
  
  # Pass a numeric vector to the cadra_search and expect that it returns an error
  expect_error(ks_genescore_mat(1:10))
  
})


test_that("ks_genescore_mat generates results consistent with ks.test() function output ", {
  

  
  res.ks.test <- suppressWarnings(ks.test(x=1:10, c(1,3,9), alternative="less"))
  mat <- matrix(c(1,0,1,0,0,0,0,0,1,0), nrow=1)
  res.cadra<- ks_genescore_mat(mat, weight=NULL, alt="less")
  
  expect_equal(res.cadra[1,1], unname(res.ks.test$statistic ))
  expect_equal(res.cadra[2,1], unname(res.ks.test$p.value ))
  
})
