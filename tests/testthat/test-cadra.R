test_that("CaDrA returns expected result ",{
  
  # Load pre-computed feature set
  data(sim_FS)
  
  # Load pre-computed input-score
  data(sim_Scores)
  
  # cadra_result <- CaDrA(
  #   FS = sim_FS, 
  #   input_score = sim_Scores, 
  #   method = "ks_pval", 
  #   weight = NULL,
  #   alternative = "less", 
  #   top_N = 1,
  #   search_start = NULL, 
  #   search_method = "both", 
  #   max_size = 7, 
  #   n_perm = 100,
  #   plot = TRUE, 
  #   smooth = TRUE, 
  #   obs_best_score = NULL,
  #   ncores = 1, 
  #   cache_path = NULL
  # )
  # 
  # testthat::expect_type(cadra_result, "list")
  # testthat::expect_length(cadra_result, 4)
  # testthat::expect_type(cadra_result[[1]], "list")
  # testthat::expect_type(cadra_result[[2]], "double")
  # testthat::expect_type(cadra_result[[3]], "double")
  # testthat::expect_type(cadra_result[[4]], "double")
  # testthat::expect_length(cadra_result[[2]], 1000)
  # testthat::expect_length(cadra_result[[3]], 1)
  # testthat::expect_length(cadra_result[[4]], 1)
  
})


