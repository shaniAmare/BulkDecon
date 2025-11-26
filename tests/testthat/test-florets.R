test_that("florets runs with basic input", {
  x <- runif(5)
  y <- runif(5)
  b <- matrix(runif(15), nrow = 3, dimnames=list(c("A","B","C"), NULL))
  expect_silent(florets(x, y, b))
})
