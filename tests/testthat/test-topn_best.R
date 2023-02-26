test_that("topn_best returns expected results",{
  
  data(topn_list)
  result <- topn_best(topn_list = topn_list)
  
  
  testthat::expect_length(result, 3L)
  testthat::expect_type(result, "list")
  testthat::expect_named(result, 
            c("feature_set", "input_score", "score"))

  testthat::expect_s4_class(result$feature_set, "SummarizedExperiment")
  testthat::expect_identical(dim(result$feature_set), c(10L, 100L))
  testthat::expect_named(result$feature_set, 
                         c("TP_8", "TP_9", "TN_216", "TP_10", "TP_6",
                           "TN_182", "TN_12", "TN_105", "TN_275", "TN_472"))
  
  
  testthat::expect_type(result$input_score, "double")
  testthat::expect_length(result$input_score, 100L)

  testthat::expect_type(result$score, "double")
  testthat::expect_length(result$score, 1L)
  testthat::expect_equal(round(result$score,5),19.50687) 
  
 })





