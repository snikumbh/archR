context("Matrix sparsification")

test_that("Handling empty matrix", {

  testMat <- matrix()
  expect_error(sparsify_mat(testMat), "empty")

  testMat <- list()
  expect_error(sparsify_mat(testMat), "Expecting matrix")
})

test_that("Handling (un)matrix", {

  testMat <- list()
  expect_error(sparsify_mat(testMat), "Expecting matrix")
})

test_that("Non-negative threshold", {
  testMat <- matrix(rnorm(20), nrow=4)
  threshold <- -0.9
  expect_error(sparsify_mat(testMat, threshold), "non-negative")

})

# TO-DO
# test_that("Sparsification itself", {
#   testMat <- matrix(rnorm(20), nrow=4)
#   threshold <- -0.9
#   expect_error(sparsify_mat(testMat, threshold), "non-negative")
#
# })
