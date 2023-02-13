test_that("prefilter_data outputs appropriate errors", {
  
  mat <- matrix(1:1000, nrow = 100)
  expect_error( prefilter_data(FS = mat) )

  
})


test_that("prefilter_data returns expected result ",{
  
  set.seed(21)
  data(sim_FS)
  # Pass a numeric vector to the cadra_search and expect that it returns an error
  result <- prefilter_data(FS = sim_FS)

  testthat::expect_s4_class(result, "SummarizedExperiment")
  expect_identical(dim(assays(result)$exprs), c(883L, 100L))

})


