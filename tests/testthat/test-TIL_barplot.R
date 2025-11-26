test_that("TIL_barplot accepts a matrix and runs", {
  mat <- matrix(runif(9), nrow = 3)
  rownames(mat) <- c("A", "B", "C")
  expect_silent(TIL_barplot(mat))
})
