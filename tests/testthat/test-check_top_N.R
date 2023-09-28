test_that("check_top_N returns expected result ",{
  
  data(sim_FS)
  data(sim_Scores)
  
  set.seed(21)
  
  rowscore <- calc_rowscore(
    FS = sim_FS,
    input_score = sim_Scores,
    meta_feature = NULL,
    method = "ks_pval",
    alternative = "less",
    weights = NULL
  )
  
  # Re-order rowscore in a decreasing order (from highest to lowest values)
  # This comes in handy when doing the top-N evaluation of
  # top N 'best' features
  rowscore <- rowscore[order(rowscore, decreasing=TRUE)]
  
  top_N_index <- check_top_N(
    rowscore = rowscore,
    top_N = 7,
    search_start = NULL,
    feature_names = rownames(sim_FS)
  )
  
  expect_type(top_N_index, "integer")
  expect_length(top_N_index, 7L)
  expect_identical(top_N_index, c(8L, 9L, 10L, 6L, 948L, 706L, 995L))
  
})


