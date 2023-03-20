test_that("ks_plot works", {


  library(SummarizedExperiment)
  data(topn_list)
  topn_best_meta <- topn_best(topn_list=topn_list)
  best_meta_feature <-  topn_best_meta[["feature_set"]]
  mat <- as.matrix(SummarizedExperiment::assay(best_meta_feature))
  or <- ifelse(colSums(mat)==0, 0, 1)
  ks_enrichment_coordinates <- ks_plot_coordinates(
     n_x = length(or),
     y = which(or==1),
     weight = NULL,
     alt = "less"
  )

  # Plot the enrichment scores
  g <- ks_plot(df = ks_enrichment_coordinates)
  expect_type(g, "list")
  expect_length(g, 9L)
  expect_s3_class(g$layers[[1]], "LayerInstance")
  expect_s3_class(g$layers[[1]]$geom, "GeomLine")
  expect_s3_class(g$layers[[2]], "LayerInstance")
  expect_s3_class(g$layers[[2]]$geom, "GeomHline")
  expect_s3_class(g$layers[[3]], "LayerInstance")
  expect_s3_class(g$layers[[3]]$geom, "GeomPoint")
  
  expect_identical(g$labels$y, "Enrichment Score (feature_set)")
})


