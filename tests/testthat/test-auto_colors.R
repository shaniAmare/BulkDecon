test_that("auto_colors returns correct length", {
  cols <- .auto_colors(5)
  expect_length(cols, 5)
})

test_that("auto_colors assigns names if provided", {
  nm <- c("A", "B", "C")
  cols <- .auto_colors(3, names = nm)
  expect_named(cols, nm)
})
