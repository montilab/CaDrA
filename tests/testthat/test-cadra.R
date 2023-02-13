test_that("CaDrA returns expected result ",{
  
  # Load pre-computed feature set
  data(sim_FS)
  
  # Load pre-computed input-score
  data(sim_Scores)
  
  set.seed(21)
  
  # Define additional parameters and start the function
  result <- CaDrA(
    FS = sim_FS, 
    input_score = sim_Scores, 
    method = "ks_pval", 
    weight = NULL, 
    alternative = "less", 
    top_N = 1,
    search_start = NULL, 
    search_method = "both", 
    max_size = 7, 
    n_perm = 10, 
    plot = FALSE, 
    smooth = TRUE, 
    obs_best_score = NULL,
    ncores = 1, 
    cache_path = NULL
  )
  
  
  testthat::expect_identical(length(result), 4L)
  testthat::expect_type(result, "list")
  testthat::expect_identical(names(result), 
              c("key","perm_best_scores","obs_best_score","perm_pval"))
  testthat::expect_type(result$key, "list")
  testthat::expect_identical(length(result$key), 12L)
  testthat::expect_identical(names(result$key), 
             c("FS", "input_score", "method", "custom_function", 
               "custom_parameters", "alternative", "weight", "top_N",
               "search_start", "search_method", "max_size", "best_score_only"))
  testthat::expect_s4_class(result$key$FS, "SummarizedExperiment")
  
  testthat::expect_identical(length(result$perm_best_scores), 10L)
  testthat::expect_equal(round(result$perm_best_scores,5), 
                         c("1.TN_84"=14.98937,
                           "2.TN_694"=16.90042,
                           "3.TN_432"=16.18984,
                           "4.TN_314"=15.41333,
                           "5.TN_636"=15.87600,
                           "6.TN_140"=14.86286,
                           "7.TN_281"=17.48639,
                           "8.TN_744"=15.51724,
                           "9.TN_504"=15.22937,
                           "10.TN_749"=14.88095))
  

})


