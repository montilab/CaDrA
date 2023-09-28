test_that("topn_best returns expected results",{
  
  data(topn_list)
  result <- topn_best(topn_list = topn_list)
  
  testthat::expect_length(result, 6L)
  testthat::expect_type(result, "list")
  testthat::expect_named(result, c("feature_set", "input_score", "score", "best_indices", "marginal_best_scores", "cumulative_best_scores"))

  testthat::expect_s4_class(result$feature_set, "SummarizedExperiment")
  testthat::expect_identical(dim(result$feature_set), c(10L, 100L))
  testthat::expect_named(result$feature_set, 
                         c('TP_7', 'TP_2', 'TP_4', 'TP_6', 'TP_3', 'TP_1', 
                           'TP_8', 'TP_9', 'TP_10', 'TP_5'))
  
  testthat::expect_type(result$input_score, "double")
  testthat::expect_length(result$input_score, 100L)

  testthat::expect_type(result$score, "double")
  testthat::expect_length(result$score, 1L)
  testthat::expect_equal(round(result$score, 5), c("TN_938"=16.66667)) 
  
})





