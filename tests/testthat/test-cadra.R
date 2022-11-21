test_that("CaDrA returns expected result ",{
  
  # Load R library
  suppressPackageStartupMessages(library(Biobase))
  
  # Load pre-computed expression set
  data(sim.ES)
  
  # set seed
  set.seed(123)
  
  # Provide a vector of continuous scores for a target profile 
  # with names to each score value 
  input_score = rnorm(n = ncol(sim.ES))
  names(input_score) <- colnames(sim.ES)
  
  cadra_result <- CaDrA(
    ES = sim.ES, input_score = input_score, 
    method = "ks", weights = NULL,
    alternative = "less", 
    metric = "pval", top_N = 1,
    search_start = NULL, 
    search_method = "both", 
    max_size = 7, n_perm = 100,
    plot = TRUE, smooth = TRUE, 
    obs_best_score = NULL,
    ncores = 1, 
    cache_path = NULL
  )
  testthat::expect_type(cadra_result, "list")
  testthat::expect_length(cadra_result, 4)
  testthat::expect_type(cadra_result[[1]], "list")
  testthat::expect_type(cadra_result[[2]], "double")
  testthat::expect_type(cadra_result[[3]], "double")
  testthat::expect_type(cadra_result[[4]], "double")
  testthat::expect_length(cadra_result[[2]], 100)
  testthat::expect_length(cadra_result[[3]], 1)
  testthat::expect_length(cadra_result[[4]], 1)
  
})


