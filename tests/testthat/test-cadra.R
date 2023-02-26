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
  
  
  testthat::expect_length(result, 4L)
  testthat::expect_type(result, "list")
  testthat::expect_named(result, 
            c("key","perm_best_scores","obs_best_score","perm_pval"))

  testthat::expect_type(result$key, "list")
  testthat::expect_length(result$key, 11L)
  testthat::expect_named(result$key, 
             c("FS", "input_score", "method", "custom_function", 
               "custom_parameters", "alternative", "weight", "top_N",
               "search_start", "search_method", "max_size"))
  testthat::expect_s4_class(result$key$FS, "SummarizedExperiment")
  
  testthat::expect_length(result$perm_best_scores, 10L)
  testthat::expect_equal(round(result$perm_best_scores,5), 
                         c("TN_84"=14.98937,
                           "TN_694"=16.90042,
                           "TN_432"=16.18984,
                           "TN_314"=15.41333,
                           "TN_636"=15.87600,
                           "TN_140"=14.86286,
                           "TN_281"=17.48639,
                           "TN_744"=15.51724,
                           "TN_504"=15.22937,
                           "TN_749"=14.88095))
 
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
  
  testthat::expect_length(result, 4L)
  testthat::expect_type(result, "list")
  testthat::expect_named(result, 
                             c("key","perm_best_scores","obs_best_score","perm_pval"))
  testthat::expect_type(result$key, "list")
  testthat::expect_length(result$key, 11L)
  testthat::expect_named(result$key, 
                             c("FS", "input_score", "method", "custom_function", 
                               "custom_parameters", "alternative", "weight", "top_N",
                               "search_start", "search_method", "max_size"))
  testthat::expect_s4_class(result$key$FS, "SummarizedExperiment")
  
  testthat::expect_length(result$perm_best_scores, 10L)
  testthat::expect_equal(round(result$perm_best_scores,2), 
                         c("TN_641"=0.97,
                           "TN_738"=0.99,
                           "TN_667"=0.99,
                           "TN_469"=0.94,
                           "TN_485"=0.96,
                           "TN_474"=0.98,
                           "TN_318"=0.98,
                           "TN_252"=0.99,
                           "TN_510"=0.95,
                           "TN_550"=0.95))
  
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
  
  
  testthat::expect_length(result, 4L)
  testthat::expect_type(result, "list")
  testthat::expect_named(result, 
            c("key","perm_best_scores","obs_best_score","perm_pval"))
  testthat::expect_type(result$key, "list")
  testthat::expect_length(result$key, 11L)
  testthat::expect_named(result$key, 
            c("FS", "input_score", "method", "custom_function", 
              "custom_parameters", "alternative", "weight", "top_N",
              "search_start", "search_method", "max_size"))
  testthat::expect_s4_class(result$key$FS, "SummarizedExperiment")
  
  testthat::expect_length(result$perm_best_scores, 10L)
  testthat::expect_equal(round(result$perm_best_scores,5), 
                         c("TN_674"=25.40974,
                           "TN_651"=29.96859,
                           "TN_490"=23.87704,
                           "TN_714"=25.97859,
                           "TN_424"=23.26229,
                           "TN_756"=30.61390,
                           "TN_845"=27.25412,
                           "TN_593"=24.18681,
                           "TN_296"=22.09454,
                           "TN_352"=23.90527))
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
  
  
  testthat::expect_length(result, 4L)
  testthat::expect_type(result, "list")
  testthat::expect_named(result, 
            c("key","perm_best_scores","obs_best_score","perm_pval"))
  testthat::expect_type(result$key, "list")
  testthat::expect_length(result$key, 11L)
  testthat::expect_named(result$key, 
            c("FS", "input_score", "method", "custom_function", 
              "custom_parameters", "alternative", "weight", "top_N",
              "search_start", "search_method", "max_size"))
  testthat::expect_s4_class(result$key$FS, "SummarizedExperiment")
  
  testthat::expect_length(result$perm_best_scores, 10L)
  testthat::expect_equal(result$perm_best_scores, 
                         c("TN_446"=2154,
                           "TN_441"=2121,
                           "TN_791"=2142,
                           "TN_691"=2132,
                           "TN_496"=2113,
                           "TN_774"=1999,
                           "TN_688"=2133,
                           "TN_247"=2158,
                           "TN_891"=2145,
                           "TN_691"=2011))
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
  
  
  testthat::expect_length(result, 4L)
  testthat::expect_type(result, "list")
  testthat::expect_named(result, 
            c("key","perm_best_scores","obs_best_score","perm_pval"))
  testthat::expect_type(result$key, "list")
  testthat::expect_length(result$key, 11L)
  testthat::expect_named(result$key, 
            c("FS", "input_score", "method", "custom_function", 
              "custom_parameters", "alternative", "weight", "top_N",
              "search_start", "search_method", "max_size"))
  testthat::expect_s4_class(result$key$FS, "SummarizedExperiment")
  
  testthat::expect_length(result$perm_best_scores, 10L)
  testthat::expect_equal(round(result$perm_best_scores,7), 
                         c("TN_607"=0.6187632,
                           "TN_651"=0.6875289,
                           "TN_490"=0.6527741,
                           "TN_482"=0.6647640,
                           "TN_424"=0.6744611,
                           "TN_888"=0.7160072,
                           "TN_845"=0.6694374,
                           "TN_128"=0.6329547,
                           "TN_432"=0.6196864,
                           "TN_282"=0.6366899))
  testthat::expect_equal(round(result$obs_best_score,5), c("TN_716"=0.68911))
  testthat::expect_equal(round(result$perm_pval,6), c(0.181818))
  
 })





