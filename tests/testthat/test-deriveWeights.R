test_that("deriveWeights runs with minimal input", {
  set.seed(123)
  norm <- matrix(runif(50), nrow = 10)
  rownames(norm) <- paste0("g", 1:10)
  colnames(norm) <- paste0("s", 1:5)

  w <- deriveWeights(norm)

  expect_true(is.matrix(w))
  expect_equal(dim(w), dim(norm))
  expect_false(any(is.na(w)))
})

test_that("deriveWeights works with raw counts", {
  set.seed(123)
  norm <- matrix(abs(rnorm(30)), nrow = 10)
  raw  <- matrix(rpois(30, 20), nrow = 10)
  rownames(norm) <- rownames(raw) <- paste0("g", 1:10)
  colnames(norm) <- colnames(raw) <- paste0("s", 1:3)

  w <- deriveWeights(norm, raw)

  expect_equal(dim(w), dim(norm))
  expect_true(all(w > 0))
})

test_that("deriveWeights aligns mismatched raw and norm", {
  set.seed(123)
  norm <- matrix(runif(40), nrow = 10)
  rownames(norm) <- paste0("gene", 1:10)
  colnames(norm) <- paste0("sample", 1:4)

  raw <- matrix(rpois(60, 20), nrow = 15)  # extra genes
  rownames(raw) <- paste0("gene", 1:15)
  colnames(raw) <- colnames(norm)

  w <- deriveWeights(norm, raw)

  expect_equal(dim(w), dim(norm))
})

test_that("deriveWeights errors when norm not matrix", {
  expect_error(deriveWeights(norm = 5), "must be a matrix")
})
