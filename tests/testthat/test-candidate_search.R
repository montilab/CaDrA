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
  ks <- candidate_search(
    FS = sim_FS, 
    input_score = sim_Scores, 
    method = "ks", 
    alternative = "less", 
    weight = NULL, 
    metric = "pval", 
    search_start = NULL, 
    search_method = "both", 
    max_size = 7, 
    best_score_only = FALSE
  )
  
  
  testthat::expect_type(ks, "list")
  testthat::expect_length(ks, 1)
  testthat::expect_type(ks[[1]], "list")
  testthat::expect_length(ks[[1]], 5)
  testthat::expect_snapshot(ks)
  
  # Run candidate_search with wilcox method
  wilcox <- candidate_search(
    FS = sim_FS, 
    input_score = sim_Scores, 
    method = "wilcox",
    alternative = "less", 
    weight = NULL, metric = "pval",
    search_start = NULL, 
    search_method = "both",
    max_size = 7, 
    best_score_only = FALSE
  )
  
  testthat::expect_type(wilcox, "list")
  testthat::expect_length(wilcox, 1)
  testthat::expect_type(wilcox[[1]], "list")
  testthat::expect_length(wilcox[[1]], 5)
  testthat::expect_snapshot(wilcox)
  
  # Run candidate_search with revealer method
  revealer <- suppressWarnings(
    candidate_search(
      FS = sim_FS, 
      input_score = sim_Scores, 
      method = "revealer",
      alternative = "less", 
      metric = "pval",
      weight = NULL,
      search_start = NULL,
      search_method = "both", 
      max_size = 7,
      best_score_only = FALSE
    )
  )
  
  testthat::expect_type(revealer, "list")
  testthat::expect_length(revealer, 1)
  testthat::expect_type(revealer[[1]], "list")
  testthat::expect_length(revealer[[1]], 5)
  testthat::expect_snapshot(revealer)
  
})


