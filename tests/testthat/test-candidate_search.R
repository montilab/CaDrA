test_that("candidate_search receives correct input arguments ", {
  
  # Pass a numeric vector to the cadra_search and expect that it returns an error
  expect_error(candidate_search(1:10))
  
})


test_that("candidate_search returns expected result ",{
  
  # Load pre-computed feature set
  data(sim_FS)
  
  # Load pre-computed input-score
  data(sim_Scores)
  
  # Define additional parameters and run the function
  result <- candidate_search(
    FS = sim_FS, 
    input_score = sim_Scores, 
    method = "ks_pval", 
    alternative = "less", 
    weight = NULL, 
    search_start = NULL, 
    search_method = "both", 
    max_size = 7, 
    best_score_only = FALSE
  )
  testthat::expect_length(result, 1L)
  testthat::expect_type(result, "list")
  testthat::expect_length(result[[1]], 3L)
  testthat::expect_s4_class(result[[1]][[1]], "SummarizedExperiment")
  
  testthat::expect_length(result[[1]][[2]], 100L)
  testthat::expect_equal(round(result[[1]][[2]][c(1:3, 100)],9),
                         c("1"=2.187332993, "2"=2.168955965, 
                           "3"=2.050084686, "100"=-2.309168876))
  
  testthat::expect_length(result[[1]][[3]], 1L)
  testthat::expect_equal(round(result[[1]][[3]],5), c("TN_716"=14.90173))
  

  # Run candidate_search with wilcox method
  result <- candidate_search(
    FS = sim_FS, 
    input_score = sim_Scores, 
    method = "wilcox_pval",
    alternative = "less", 
    weight = NULL, 
    search_start = NULL, 
    search_method = "both",
    max_size = 7, 
    best_score_only = FALSE
  )
  
  testthat::expect_length(result, 1L)
  testthat::expect_type(result, "list")
  testthat::expect_length(result[[1]], 3L)
  testthat::expect_s4_class(result[[1]][[1]], "SummarizedExperiment")
  
  testthat::expect_length(result[[1]][[2]], 100L)
  testthat::expect_equal(round(result[[1]][[2]][c(1:3, 100)],6),
                         c("1"=2.187333, "2"=2.168956, 
                           "3"=2.050085, "100"=-2.309169))
  
  testthat::expect_length(result[[1]][[3]], 1L)
  testthat::expect_equal(round(result[[1]][[3]],5), c("TP_9"=27.75113))
  
  
  
  # Run candidate_search with revealer method
  result <- suppressWarnings(
    candidate_search(
      FS = sim_FS, 
      input_score = sim_Scores, 
      method = "revealer",
      alternative = "less", 
      weight = NULL,
      search_start = NULL,
      search_method = "both", 
      max_size = 7,
      best_score_only = FALSE
    )
  )
  testthat::expect_length(result, 1L)
  testthat::expect_type(result, "list")
  testthat::expect_length(result[[1]], 3L)
  testthat::expect_s4_class(result[[1]][[1]], "SummarizedExperiment")
  
  testthat::expect_length(result[[1]][[2]], 100L)
  testthat::expect_equal(round(result[[1]][[2]][c(1:3, 100)],6),
                         c("1"=2.187333, "2"=2.168956, 
                           "3"=2.050085, "100"=-2.309169))
  
  testthat::expect_length(result[[1]][[3]], 1L)
  testthat::expect_equal(round(result[[1]][[3]],5), c("TN_716"=0.68911))
  

  
})


