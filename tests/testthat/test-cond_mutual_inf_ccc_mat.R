test_that("cond_mutual_inf_ccc_mat returns expected result", {

  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
  ZM <- cbind(mutual_info_df$Yc, mutual_info_df$Wc)
  result <- cond_mutual_inf_ccc_mat(mutual_info_df$Zc_XcYcWc, M, ZM)

  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.1171533, 0.2192397), tolerance = 0.000001)

})

test_that("cond_mutual_inf_ccc_mat issues error messages when vector and matrix have different sizes", {

  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
  ZM <- cbind(mutual_info_df$Yc, mutual_info_df$Wc)
  expect_error( cond_mutual_inf_ccc_mat(mutual_info_df$Zc_XcYcWc[-1], M, ZM))
  expect_error( cond_mutual_inf_ccc_mat(mutual_info_df$Zc_XcYcWc,
                                       M, ZM[-1,]))
  expect_error( cond_mutual_inf_ccc_mat(mutual_info_df$Zc_XcYcWc,
                                       M, ZM[, 1]))
  expect_error( cond_mutual_inf_ccc_mat(mutual_info_df$Zc_XcYcWc,
                                       M[-1,], ZM))
  expect_error( cond_mutual_inf_ccc_mat(mutual_info_df$Zc_XcYcWc,
                                       M[, 1], ZM))
})


test_that("cond_mutual_inf_ccc_mat issues error messages when the value of k is too large", {

  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
  ZM <- cbind(mutual_info_df$Yc, mutual_info_df$Wc)
  expect_error( cond_mutual_inf_ccc_mat(mutual_info_df$Zc_XcYcWc, M, ZM, k=101))
})

