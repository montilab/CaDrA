test_that("meta_plot works", {
  # Load pre-computed Top-N list generated for sim_FS dataset
  data(topn_list)

  # With the results obtained from top-N evaluation,
  # We can find the combination of features that gives the best score in
  # top N searches
  topn_best_meta <- topn_best(topn_list = topn_list)

  # Now we can plot this set of best meta-feature
  expect_silent(  meta_plot(topn_best_list = topn_best_meta))
  
})
