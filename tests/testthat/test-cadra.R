test_that("CaDrA returns expected result for ks algorithm",{
  
  # Load pre-computed feature set
  data(sim_FS)
  
  # Load pre-computed input-score
  data(sim_Scores)
  
  # Set seed
  set.seed(21)
  
  # ks_pval
  result <- CaDrA(
    FS = sim_FS, 
    input_score = sim_Scores, 
    method = "ks_pval", 
    weights = NULL, 
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
    cache = FALSE,
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
               "custom_parameters", "alternative", "weights", "top_N",
               "search_start", "search_method", "max_size"))
  testthat::expect_type(result$key$FS, "double")
  
  testthat::expect_equal(round(result$perm_best_scores[1:10], 5), 
                         c('TN_780'=3.85692,
                           'TN_989'=4.24386,
                           'TN_974'=4.64671,
                           'TN_271'=4.46158,
                           'TN_427'=5.7888,
                           'TN_983'=4.83556,
                           'TN_170'=4.59927,
                           'TN_669'=4.40737,
                           'TN_605'=4.50667,
                           'TN_874'=3.60332))
 
  testthat::expect_equal(round(result$obs_best_score, 5), c("TP_8"=14.13128))
  
  # A smooth factor of 1
  c <- 1
  
  # Add a smoothing factor of 1 
  # This is just to not return a p-value of 0
  testthat::expect_equal(
    round((sum(result$perm_best_scores[1:10] > result$obs_best_score)+c)/(10+c), 7), 
    c(0.0909091)
  )
  
  # # Set seed
  # set.seed(21)
  # 
  # # ks_score
  # result <- CaDrA(
  #   FS = sim_FS,
  #   input_score = sim_Scores,
  #   method = "ks_score",
  #   weights = NULL,
  #   alternative = "less",
  #   top_N = 1,
  #   search_start = NULL,
  #   search_method = "both",
  #   max_size = 7,
  #   n_perm = 10,
  #   plot = FALSE,
  #   smooth = TRUE,
  #   obs_best_score = NULL,
  #   ncores = 1,
  #   cache = FALSE,
  #   cache_path = NULL
  # )
  # 
  # testthat::expect_length(result, 4L)
  # testthat::expect_type(result, "list")
  # testthat::expect_named(result,
  #                            c("key","perm_best_scores","obs_best_score","perm_pval"))
  # testthat::expect_type(result$key, "list")
  # testthat::expect_length(result$key, 11L)
  # testthat::expect_named(result$key,
  #                            c("FS", "input_score", "method", "custom_function",
  #                              "custom_parameters", "alternative", "weights", "top_N",
  #                              "search_start", "search_method", "max_size"))
  # testthat::expect_type(result$key$FS, "double")
  # 
  # testthat::expect_equal(round(result$perm_best_scores[1:10], 5),
  #                        c('TP_8'=0.34,
  #                          'TP_10'=0.54,
  #                          'TP_9'=0.37,
  #                          'TP_6'=0.40,
  #                          'TP_9'=0.38,
  #                          'TP_9'=0.52,
  #                          'TP_2'=0.44,
  #                          'TP_2'=0.40,
  #                          'TP_4'=0.44,
  #                          'TP_9'=0.49))
  # 
  # testthat::expect_equal(round(result$obs_best_score, 2), c("TP_9"=0.66))
  # 
  # # A smooth factor of 1
  # c <- 1
  # 
  # # Add a smoothing factor of 1
  # # This is just to not return a p-value of 0
  # testthat::expect_equal(
  #   round((sum(result$perm_best_scores[1:10] > result$obs_best_score)+c)/(10+c), 6),
  #   c(0.090909)
  # )
  
})

# ========================================================================= #
test_that("CaDrA returns expected result for Wilcoxon algorithm",{
  
  # Load pre-computed feature set
  data(sim_FS)
  
  # Load pre-computed input-score
  data(sim_Scores)
  
  # Set seed
  set.seed(21)
  
  # wilcox_pval
  result <- CaDrA(
    FS = sim_FS, 
    input_score = sim_Scores, 
    method = "wilcox_pval", 
    weights = NULL, 
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
    cache = FALSE,
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
              "custom_parameters", "alternative", "weights", "top_N",
              "search_start", "search_method", "max_size"))
  testthat::expect_type(result$key$FS, "double")
  
  testthat::expect_equal(round(result$perm_best_scores[1:10], 5), 
                         c('TN_780'=14.98171,
                           'TN_97'=16.91855,
                           'TN_974'=22.01412,
                           'TN_271'=19.23964,
                           'TN_141'=19.15244,
                           'TN_315'=19.56755,
                           'TN_689'=16.51065,
                           'TN_413'=20.50307,
                           'TN_940'=20.42022,
                           'TN_894'=17.80290))
  testthat::expect_equal(round(result$obs_best_score, 5), c("TN_129"=21.35299))
  
  # A smooth factor of 1
  c <- 1
  
  # Add a smoothing factor of 1 
  # This is just to not return a p-value of 0
  testthat::expect_equal(
    round((sum(result$perm_best_scores[1:10] > result$obs_best_score)+c)/(10+c), 6), 
    c(0.181818)
  )
  
  # Set seed
  set.seed(21)
  
  # wilcox_score
  result <- CaDrA(
    FS = sim_FS, 
    input_score = sim_Scores, 
    method = "wilcox_score", 
    weights = NULL, 
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
    cache = FALSE,
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
              "custom_parameters", "alternative", "weights", "top_N",
              "search_start", "search_method", "max_size"))
  testthat::expect_type(result$key$FS, "double")
  
  testthat::expect_equal(result$perm_best_scores[1:10], 
                         c('TN_200'=1879,
                           'TN_927'=1860,
                           'TN_596'=1918,
                           'TN_977'=1875,
                           'TN_532'=1789,
                           'TN_951'=1826,
                           'TN_788'=1748,
                           'TN_884'=1904,
                           'TN_181'=1975,
                           'TN_987'=1974))
  testthat::expect_equal(result$obs_best_score, c("TN_544"=1738))
  
  # A smooth factor of 1
  c <- 1
  
  # Add a smoothing factor of 1 
  # This is just to not return a p-value of 0
  testthat::expect_equal(
    round((sum(result$perm_best_scores[1:10] > result$obs_best_score)+c)/(10+c), 6), 
    c(1)
  )
  
})


# ========================================================================= #
test_that("CaDrA returns expected result for Revealer algorithm", {
  
  # Load pre-computed feature set
  data(sim_FS)
  
  # Load pre-computed input-score
  data(sim_Scores)
  
  # Set seed
  set.seed(21)
  
  # Revealer
  result <- CaDrA(
    FS = sim_FS, 
    input_score = sim_Scores, 
    method = "revealer", 
    weights = NULL, 
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
    cache = FALSE,
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
              "custom_parameters", "alternative", "weights", "top_N",
              "search_start", "search_method", "max_size"))
  testthat::expect_type(result$key$FS, "double")
  
  testthat::expect_equal(round(result$perm_best_scores[1:10], 7), 
                         c('TN_780'=0.3901349,
                           'TN_97'=0.4065603,
                           'TN_974'=0.3892586,
                           'TN_393'=0.3805595,
                           'TN_369'=0.4212498,
                           'TN_983'=0.4050481,
                           'TN_681'=0.3577511,
                           'TN_125'=0.3747675,
                           'TN_940'=0.3845112,
                           'TN_351'=0.4654268))
  testthat::expect_equal(round(result$obs_best_score, 5), c("TN_985"=0.37856))
  
  # A smooth factor of 1
  c <- 1
  
  # Add a smoothing factor of 1 
  # This is just to not return a p-value of 0
  testthat::expect_equal(
    round((sum(result$perm_best_scores[1:10] > result$obs_best_score)+c)/(10+c), 6), 
    c(0.818182)
  )
  
})





