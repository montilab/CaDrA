test_that("mutual_inf_cc_mat returns expected result", {

  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
  result <- mutual_inf_cc_mat(mutual_info_df$Zc_XcYcWc, M)

  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.000000, 0.199844), tolerance = 0.000001)

})

test_that("mutual_inf_cc_mat issues error messages when vector and matrix have different sizes", {

  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)

  expect_error( mutual_inf_cc_mat(mutual_info_df$Zc_XcYcWc[-1], M))
})


test_that("mutual_inf_cc_mat issues error messages when the value of k is too large", {

  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)

  expect_error(  mutual_inf_cc_mat(mutual_info_df$Zc_XcYcWc, M, k=150) )
})

