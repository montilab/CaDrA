test_that("cond_mutual_inf_cdd_vec returns expected result", {

  data(mutual_info_df)

  result <- cond_mutual_inf_cdd_vec(mutual_info_df$Zc_XdYd, mutual_info_df$Xd,
                                   mutual_info_df$Yd)

  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(result, 0.1338664, tolerance = 0.00001)
})


test_that("cond_mutual_inf_cdd_vec issues error messages when vector sizes are different", {

  data(mutual_info_df)

  expect_error( cond_mutual_inf_cdd_vec(mutual_info_df$Zc_XdYd[-1],
                                       mutual_info_df$Xd,
                                       mutual_info_df$Yd))
  expect_error( cond_mutual_inf_cdd_vec(mutual_info_df$Zc_XdYd,
                                       mutual_info_df$Xd[-1],
                                       mutual_info_df$Yd))
  expect_error( cond_mutual_inf_cdd_vec(mutual_info_df$Zc_XdYd,
                                       mutual_info_df$Xd,
                                       mutual_info_df$Yd[-1]))
})


test_that("cond_mutual_inf_cdd_vec issues error messages when the value of k is too large", {

  data(mutual_info_df)

  expect_error( cond_mutual_inf_cdd_vec(mutual_info_df$Zc_XdYd,
                                       mutual_info_df$Xd,
                                       mutual_info_df$Yd, k=101))
})
