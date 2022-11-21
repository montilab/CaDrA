test_that("candidate_search receives correct input arguments ", {
  
  # Pass a numeric vector to the cadra_search and expect that it returns an error
  expect_error(candidate_search(1:10))
  
})


test_that("candidate_search returns expected result ",{
  
  # Load R library
  suppressPackageStartupMessages(library(Biobase))
  
  # Load pre-computed expression set
  data(sim.ES)
  
  # set seed
  set.seed(123)
  
  # Provide a vector of continuous scores for a target profile with 
  # names to each score value 
  input_score = rnorm(n = ncol(sim.ES))
  names(input_score) <- colnames(sim.ES)
  
  # Define additional parameters and run the function
  ks <- candidate_search(
    ES = sim.ES, 
    input_score = input_score, 
    method = "ks", 
    alternative = "less", 
    weights = NULL, 
    metric = "pval", 
    search_start = NULL, 
    search_method = "both", 
    max_size = 7, 
    best_score_only = FALSE
  )
  
  
  testthat::expect_type(ks, "list")
  testthat::expect_length(ks, 1)
  testthat::expect_type(ks[[1]], "list")
  testthat::expect_length(ks[[1]], 3)
  testthat::expect_snapshot(ks)
  
  # Run candidate_search with wilcox method
  wilcox <- candidate_search(
    ES = sim.ES, input_score = input_score, method = "wilcox",
    alternative = "less", weights = NULL, metric = "pval",
    search_start = NULL, search_method = "both",
    max_size = 7, best_score_only = FALSE
  )
  
  testthat::expect_type(wilcox, "list")
  testthat::expect_length(wilcox, 1)
  testthat::expect_type(wilcox[[1]], "list")
  testthat::expect_length(wilcox[[1]], 3)
  testthat::expect_snapshot(wilcox)
  
  # Run candidate_search with revealer method
  revealer <- suppressWarnings(candidate_search(
    ES = sim.ES, input_score = input_score, method = "revealer",
    alternative = "less", 
    metric = "pval",
    weights = NULL,
    search_start = NULL,
    search_method = "both", 
    max_size = 7,
    best_score_only = FALSE,
    verbose = FALSE
  ))
  
  testthat::expect_type(revealer, "list")
  testthat::expect_length(revealer, 1)
  testthat::expect_type(revealer[[1]], "list")
  testthat::expect_length(revealer[[1]], 3)
  testthat::expect_snapshot(revealer)
  
})


