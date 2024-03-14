test_that("generate_permutations returns expected results", {
  
  data(sim_Scores)
  set.seed(123)
  n_perm = 1000

  result <- generate_permutations(
    input_score = sim_Scores,
    n_perm = n_perm
  )
  
  testthat::expect_identical(dim(result), c(1000L, 100L))
  testthat::expect_type(result, "double")
})
