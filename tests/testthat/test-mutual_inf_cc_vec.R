test_that("mutual_inf_cc_vec returns expected result", {

  data(mutual_info_df)

  result <- mutual_inf_cc_vec(mutual_info_df$Xc, mutual_info_df$Zc_XcYc)
  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(result, 0.0, tolerance = 0.00001)


  result <- mutual_inf_cc_vec(mutual_info_df$Yc, mutual_info_df$Zc_XcYc)
  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(result, 0.2738658, tolerance = 0.00001)
})

test_that("mutual_inf_cc_vec returns expected result with seed", {

  x <- c(0.91 , 0.08, 0.13, 0.01, 0.08,
         0.25, 0.86, 0.13 , 0.26, 0.7, 0.18, 0.7, 0.26, 0.27 , 0.46)

  y <- c(1.75,  1.14,  0.99,  0.96,  1.08,
         1.18,  1.63,  1.03,  0.95, -0.31, 0.89,  1.45,  1.02,  0.97,  1.25)

  result <- mutual_inf_cc_vec(x, y, k=3, seed=654321)
  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(result, 0.3558762, tolerance = 0.00001)
})


test_that("mutual_inf_cc_vec issues error messages when vector sizes are different", {

  x <- c(0.91 , 0.08, 0.13, 0.01, 0.08,
         0.25, 0.86, 0.13 , 0.26, 0.7, 0.18, 0.7, 0.26, 0.27 , 0.46)

  y <- c(1.75,  1.14,  0.99,  0.96,  1.08)


  expect_error(  mutual_inf_cc_vec(x, y, k=3) )
})

test_that("mutual_inf_cc_vec issues error messages when the value of k is too large", {

  x <- c(0.91 , 0.08, 0.13, 0.01, 0.08,
         0.25, 0.86, 0.13 , 0.26, 0.7, 0.18, 0.7, 0.26, 0.27 , 0.46)

  y <- c(1.75,  1.14,  0.99,  0.96,  1.08,
         1.18,  1.63,  1.03,  0.95, -0.31, 0.89,  1.45,  1.02,  0.97,  1.25)


  expect_error(  mutual_inf_cc_vec(x, y, k=16) )
})
