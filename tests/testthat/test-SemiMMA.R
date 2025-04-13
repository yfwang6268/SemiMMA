test_that("SemiMMA handles basic input", {
  data(data_matrix, package = "SemiMMA")
  # Add column names
  n = data_matrix[,ncol(data_matrix)]
  y = data_matrix[,seq(1,ncol(data_matrix)-1,2)]
  s = data_matrix[,seq(2,ncol(data_matrix)-1,2)]
  s_imputed <- impute_within_study_sd(s, n)
  result <- SemiMMA(y,s,n)
  expect_equal(sum(is.na(s_imputed)), 0)
  expect_equal(length(result), 10)
})
