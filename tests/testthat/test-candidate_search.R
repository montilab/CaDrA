test_that("cadra_search receives correct input arguments ", {

  # Pass a numeric vector to the cadra_search and expect that it returns an error
  expect_error(cadra_search(1:10))
  
})