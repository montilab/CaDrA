test_that("ks_genescore_mat receives correct input arguments ", {
  
  # Pass a numeric vector to the cadra_search and expect that it returns an error
  expect_error(ks_genescore_mat(1:10))
  
})
