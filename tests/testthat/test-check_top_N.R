test_that("check_top_N returns expected result ",{
  
  library(SummarizedExperiment)
  data(sim_FS)
  data(sim_Scores)
  
  set.seed(21)
  
  rowscore <- calc_rowscore(
    FS_mat = SummarizedExperiment::assay(sim_FS),
    input_score = sim_Scores,
    method = "ks_pval",
    alternative = "less",
    weight = NULL)
    
    top_N_index <- check_top_N(
    rowscore = rowscore,
    top_N = 7,
    search_start = NULL,
    feature_names = rownames(sim_FS))
    

    expect_type(top_N_index, "integer")
    expect_length(top_N_index, 7L)
    expect_identical(top_N_index, c(726L, 8L, 9L, 854L, 973L, 589L, 118L))
    
})


