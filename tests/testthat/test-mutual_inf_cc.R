test_that("mutual_inf_cc returns expected result", {
  set.seed(654321)
  
  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
  tM <- t(M)
  result <- mutual_inf_cc(mutual_info_df$Zc_XcYcWc, tM)
  
  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.000000, 0.199844), tolerance = 0.000001)
  
  result <- mutual_inf_cc(mutual_info_df$Zc_XcYcWc, t(M), k=4)
  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.01247941, 0.14108933), tolerance = 0.000001)  
  
})

test_that("mutual_inf_cc issues error messages when vector and matrix have different sizes", {
  
  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
  
  expect_error( mutual_inf_cc(mutual_info_df$Zc_XcYcWc[-1], t(M)))
})


test_that("mutual_inf_cc issues error messages when the value of k is too large", {
  
  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
  
  expect_error(  mutual_inf_cc(mutual_info_df$Zc_XcYcWc, t(M), k=150) )
})

test_that("mutual_inf_cc returns expected result", {
  set.seed(654321)
  
  data(mutual_info_df)
  
  result <- mutual_inf_cc(mutual_info_df$Xc, t(mutual_info_df$Zc_XcYc))
  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(result, 0.0, tolerance = 0.00001)
  
  result <- mutual_inf_cc(mutual_info_df$Xc, t(mutual_info_df$Zc_XcYc),k=5)
  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(result, 0.010756, tolerance = 0.00001)
  
  result <- mutual_inf_cc(mutual_info_df$Yc, t(mutual_info_df$Zc_XcYc))
  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(result, 0.2738658, tolerance = 0.00001)
})

test_that("mutual_inf_cc returns expected result with seed", {
  
  x <- c(0.91 , 0.08, 0.13, 0.01, 0.08,
         0.25, 0.86, 0.13 , 0.26, 0.7, 0.18, 0.7, 0.26, 0.27 , 0.46)
  
  y <- c(1.75,  1.14,  0.99,  0.96,  1.08,
         1.18,  1.63,  1.03,  0.95, -0.31, 0.89,  1.45,  1.02,  0.97,  1.25)
  
  set.seed(654321)
  result <- mutual_inf_cc(x, y, k=3)
  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(result, 0.3558762, tolerance = 0.00001)
})


test_that("mutual_inf_cc issues error messages when vector sizes are different", {
  
  x <- c(0.91 , 0.08, 0.13, 0.01, 0.08,
         0.25, 0.86, 0.13 , 0.26, 0.7, 0.18, 0.7, 0.26, 0.27 , 0.46)
  
  y <- c(1.75,  1.14,  0.99,  0.96,  1.08)
  
  
  expect_error(  mutual_inf_cc(x, y, k=3) )
})

 
