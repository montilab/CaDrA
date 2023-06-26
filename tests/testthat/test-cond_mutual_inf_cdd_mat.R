test_that("cond_mutual_inf_cdd_mat returns expected result", {

  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  ZM <- cbind(mutual_info_df$Yd, mutual_info_df$Wd)
  result <- cond_mutual_inf_cdd_mat(mutual_info_df$Zc_XdYdWd, M=M, Z=ZM)

  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.1757598, 0.1086227), tolerance = 0.000001)

})

test_that("cond_mutual_inf_cdd_mat issues error messages when vector and matrix have different sizes", {

  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  ZM <- cbind(mutual_info_df$Yd, mutual_info_df$Wd)
  expect_error( cond_mutual_inf_cdd_mat(mutual_info_df$Zc_XdYdWd[-1],
                                       M,
                                       ZM))
  expect_error( cond_mutual_inf_ccc_mat(mutual_info_df$Zc_XdYdWd,
                                       M, ZM[-1,]))
  expect_error( cond_mutual_inf_ccc_mat(mutual_info_df$Zc_XdYdWd,
                                       M, ZM[, 1]))
  expect_error( cond_mutual_inf_ccc_mat(mutual_info_df$Zc_XdYdWd,
                                       M[-1,], ZM))
  expect_error( cond_mutual_inf_ccc_mat(mutual_info_df$Zc_XdYdWd,
                                       M[, 1], ZM))
})


test_that("cond_mutual_inf_cdd_mat issues error messages when the value of k is too large", {

  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  ZM <- cbind(mutual_info_df$Yd, mutual_info_df$Wd)
  expect_error( cond_mutual_inf_cdd_mat(mutual_info_df$Zc_XdYdWd[-1],
                                       M,
                                       ZM, k=101))
})

