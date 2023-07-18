test_that("topn_plot works", {

  data(topn_list)
  
  # Get the overlapping top N plot
  expect_silent(topn_plot(topn_list = topn_list))

})


