test_that("CaDrA returns expected result for ks algorithm",{
  
  # Load pre-computed feature set
  data(sim_FS)
  
  # Load pre-computed input-score
  data(sim_Scores)
  
  set.seed(21)
  
  # ks_pval
  result <- CaDrA(
    FS = sim_FS, 
    input_score = sim_Scores, 
    method = "ks_pval", 
    weight = NULL, 
    alternative = "less", 
    top_N = 1,
    search_start = NULL, 
    search_method = "both", 
    max_size = 7, 
    n_perm = 10, 
    plot = FALSE, 
    smooth = TRUE, 
    obs_best_score = NULL,
    ncores = 1, 
    cache_path = NULL
  )
  
  
  testthat::expect_identical(length(result), 4L)
  testthat::expect_type(result, "list")
  testthat::expect_identical(names(result), 
              c("key","perm_best_scores","obs_best_score","perm_pval"))
  testthat::expect_type(result$key, "list")
  testthat::expect_identical(length(result$key), 12L)
  testthat::expect_identical(names(result$key), 
             c("FS", "input_score", "method", "custom_function", 
               "custom_parameters", "alternative", "weight", "top_N",
               "search_start", "search_method", "max_size", "best_score_only"))
  testthat::expect_s4_class(result$key$FS, "SummarizedExperiment")
  
  testthat::expect_identical(length(result$perm_best_scores), 10L)
  testthat::expect_equal(round(result$perm_best_scores,5), 
                         c("1.TN_84"=14.98937,
                           "2.TN_694"=16.90042,
                           "3.TN_432"=16.18984,
                           "4.TN_314"=15.41333,
                           "5.TN_636"=15.87600,
                           "6.TN_140"=14.86286,
                           "7.TN_281"=17.48639,
                           "8.TN_744"=15.51724,
                           "9.TN_504"=15.22937,
                           "10.TN_749"=14.88095))
 
  testthat::expect_equal(round(result$obs_best_score,5), c("TN_716"=14.90173))
  testthat::expect_equal(round(result$perm_pval,7), c(0.8181818))
  
   
  set.seed(21)
  # ks_score
  result <- CaDrA(
    FS = sim_FS, 
    input_score = sim_Scores, 
    method = "ks_score", 
    weight = NULL, 
    alternative = "less", 
    top_N = 1,
    search_start = NULL, 
    search_method = "both", 
    max_size = 7, 
    n_perm = 10, 
    plot = FALSE, 
    smooth = TRUE, 
    obs_best_score = NULL,
    ncores = 1, 
    cache_path = NULL
  )
  
  testthat::expect_identical(length(result), 4L)
  testthat::expect_type(result, "list")
  testthat::expect_identical(names(result), 
                             c("key","perm_best_scores","obs_best_score","perm_pval"))
  testthat::expect_type(result$key, "list")
  testthat::expect_identical(length(result$key), 12L)
  testthat::expect_identical(names(result$key), 
                             c("FS", "input_score", "method", "custom_function", 
                               "custom_parameters", "alternative", "weight", "top_N",
                               "search_start", "search_method", "max_size", "best_score_only"))
  testthat::expect_s4_class(result$key$FS, "SummarizedExperiment")
  
  testthat::expect_identical(length(result$perm_best_scores), 10L)
  testthat::expect_equal(round(result$perm_best_scores,2), 
                         c("1.TN_641"=0.97,
                           "2.TN_738"=0.99,
                           "3.TN_667"=0.99,
                           "4.TN_469"=0.94,
                           "5.TN_485"=0.96,
                           "6.TN_474"=0.98,
                           "7.TN_318"=0.98,
                           "8.TN_252"=0.99,
                           "9.TN_510"=0.95,
                           "10.TN_550"=0.95))
  
  testthat::expect_equal(round(result$obs_best_score,2), c("TN_278"=0.98))
  testthat::expect_equal(round(result$perm_pval,6), c(0.363636))
  
})

# ========================================================================= #
test_that("CaDrA returns expected result for Wilcoxon algorithm",{
  
  # Load pre-computed feature set
  data(sim_FS)
  
  # Load pre-computed input-score
  data(sim_Scores)
  
  set.seed(21)
  
  # wilcox_pval
  result <- CaDrA(
    FS = sim_FS, 
    input_score = sim_Scores, 
    method = "wilcox_pval", 
    weight = NULL, 
    alternative = "less", 
    top_N = 1,
    search_start = NULL, 
    search_method = "both", 
    max_size = 7, 
    n_perm = 10, 
    plot = FALSE, 
    smooth = TRUE, 
    obs_best_score = NULL,
    ncores = 1, 
    cache_path = NULL
  )
  
  
  testthat::expect_identical(length(result), 4L)
  testthat::expect_type(result, "list")
  testthat::expect_identical(names(result), 
                             c("key","perm_best_scores","obs_best_score","perm_pval"))
  testthat::expect_type(result$key, "list")
  testthat::expect_identical(length(result$key), 12L)
  testthat::expect_identical(names(result$key), 
                             c("FS", "input_score", "method", "custom_function", 
                               "custom_parameters", "alternative", "weight", "top_N",
                               "search_start", "search_method", "max_size", "best_score_only"))
  testthat::expect_s4_class(result$key$FS, "SummarizedExperiment")
  
  testthat::expect_identical(length(result$perm_best_scores), 10L)
  testthat::expect_equal(round(result$perm_best_scores,5), 
                         c("1.TN_674"=25.40974,
                           "2.TN_651"=29.96859,
                           "3.TN_490"=23.87704,
                           "4.TN_714"=25.97859,
                           "5.TN_424"=23.26229,
                           "6.TN_756"=30.61390,
                           "7.TN_845"=27.25412,
                           "8.TN_593"=24.18681,
                           "9.TN_296"=22.09454,
                           "10.TN_352"=23.90527))
  testthat::expect_equal(round(result$obs_best_score,5), c("TP_9"=27.75113))
  testthat::expect_equal(round(result$perm_pval,6), c(0.272727))
  
  
  set.seed(21)
  
  # wilcox_score
  result <- CaDrA(
    FS = sim_FS, 
    input_score = sim_Scores, 
    method = "wilcox_score", 
    weight = NULL, 
    alternative = "less", 
    top_N = 1,
    search_start = NULL, 
    search_method = "both", 
    max_size = 7, 
    n_perm = 10, 
    plot = FALSE, 
    smooth = TRUE, 
    obs_best_score = NULL,
    ncores = 1, 
    cache_path = NULL
  )
  
  
  testthat::expect_identical(length(result), 4L)
  testthat::expect_type(result, "list")
  testthat::expect_identical(names(result), 
                             c("key","perm_best_scores","obs_best_score","perm_pval"))
  testthat::expect_type(result$key, "list")
  testthat::expect_identical(length(result$key), 12L)
  testthat::expect_identical(names(result$key), 
                             c("FS", "input_score", "method", "custom_function", 
                               "custom_parameters", "alternative", "weight", "top_N",
                               "search_start", "search_method", "max_size", "best_score_only"))
  testthat::expect_s4_class(result$key$FS, "SummarizedExperiment")
  
  testthat::expect_identical(length(result$perm_best_scores), 10L)
  testthat::expect_equal(result$perm_best_scores, 
                         c("1.TN_446"=2154,
                           "2.TN_441"=2121,
                           "3.TN_791"=2142,
                           "4.TN_691"=2132,
                           "5.TN_496"=2113,
                           "6.TN_774"=1999,
                           "7.TN_688"=2133,
                           "8.TN_247"=2158,
                           "9.TN_891"=2145,
                           "10.TN_691"=2011))
  testthat::expect_equal(result$obs_best_score, c("TN_277"=2150))
  testthat::expect_equal(round(result$perm_pval,6), c(0.272727))
})


# ========================================================================= #
test_that("CaDrA returns expected result for Revealer algorithm",{
  
  # Load pre-computed feature set
  data(sim_FS)
  
  # Load pre-computed input-score
  data(sim_Scores)
  
  set.seed(21)
  
  # Revealer
  result <- CaDrA(
    FS = sim_FS, 
    input_score = sim_Scores, 
    method = "revealer", 
    weight = NULL, 
    alternative = "less", 
    top_N = 1,
    search_start = NULL, 
    search_method = "both", 
    max_size = 7, 
    n_perm = 10, 
    plot = FALSE, 
    smooth = TRUE, 
    obs_best_score = NULL,
    ncores = 1, 
    cache_path = NULL
  )
  
  
  testthat::expect_identical(length(result), 4L)
  testthat::expect_type(result, "list")
  testthat::expect_identical(names(result), 
                             c("key","perm_best_scores","obs_best_score","perm_pval"))
  testthat::expect_type(result$key, "list")
  testthat::expect_identical(length(result$key), 12L)
  testthat::expect_identical(names(result$key), 
                             c("FS", "input_score", "method", "custom_function", 
                               "custom_parameters", "alternative", "weight", "top_N",
                               "search_start", "search_method", "max_size", "best_score_only"))
  testthat::expect_s4_class(result$key$FS, "SummarizedExperiment")
  
  testthat::expect_identical(length(result$perm_best_scores), 10L)
  testthat::expect_equal(round(result$perm_best_scores,7), 
                         c("1.TN_607"=0.6187632,
                           "2.TN_651"=0.6875289,
                           "3.TN_490"=0.6527741,
                           "4.TN_482"=0.6647640,
                           "5.TN_424"=0.6744611,
                           "6.TN_888"=0.7160072,
                           "7.TN_845"=0.6694374,
                           "8.TN_128"=0.6329547,
                           "9.TN_432"=0.6196864,
                           "10.TN_282"=0.6366899))
  testthat::expect_equal(round(result$obs_best_score,5), c("TN_716"=0.68911))
  testthat::expect_equal(round(result$perm_pval,6), c(0.181818))
  
 })





