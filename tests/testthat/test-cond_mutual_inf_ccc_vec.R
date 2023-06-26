test_that("cond_mutual_inf_ccc_vec returns expected result", {

  data(mutual_info_df)

  result <- cond_mutual_inf_ccc_vec(mutual_info_df$Zc_XcYc,
                                   mutual_info_df$Xc, mutual_info_df$Yc)
  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(result, 0.2936858, tolerance = 0.00001)

})


test_that("cond_mutual_inf_ccc_vec issues error messages when vector sizes are different", {

  data(mutual_info_df)

  expect_error(   cond_mutual_inf_ccc_vec(mutual_info_df$Zc_XcYc[-1],
                                   mutual_info_df$Xc, mutual_info_df$Yc) )
  expect_error(   cond_mutual_inf_ccc_vec(mutual_info_df$Zc_XcYc,
                                         mutual_info_df$Xc[-1], mutual_info_df$Yc) )
  expect_error(   cond_mutual_inf_ccc_vec(mutual_info_df$Zc_XcYc,
                                         mutual_info_df$Xc, mutual_info_df$Yc[-1]) )

})

test_that("cond_mutual_inf_ccc_vec issues error messages when the value of k is too large", {

  expect_error(   cond_mutual_inf_ccc_vec(mutual_info_df$Zc_XcYc,
                                         mutual_info_df$Xc,
                                         mutual_info_df$Yc, k=101) )
})
