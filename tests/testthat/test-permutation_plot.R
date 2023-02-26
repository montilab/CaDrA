test_that("permutation_plot works", {


  data(perm_res)

  # Plot the permutation results
  g <- permutation_plot(perm_res)
  
  expect_type(g, "list")
  expect_length(g, 9L)
  expect_s3_class(g$layers[[1]], "LayerInstance")
  expect_s3_class(g$layers[[1]]$geom, "GeomBar")
  expect_s3_class(g$layers[[2]], "LayerInstance")
  expect_s3_class(g$layers[[2]]$geom, "GeomVline")
  
  expect_identical(g$labels$x, "Score")
  expect_identical(g$labels$y, "Count")
})


