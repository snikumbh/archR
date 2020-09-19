context("One-hot encoding functions")

test_that("Report empty sequence", {
  #
  testSeqsNotOK <- c()
  expect_error(get_one_hot_encoded_seqs(testSeqsNotOK), "Empty or NULL found")
  expect_error(.one_hot_encode_sinuc(testSeqsNotOK), "Empty or NULL found")
})

test_that("One-hot encoding is correctly done", {
  testSeqs <- c("GGCT", "ATCG")
  testSeqsB <- list(c("G", "G", "C", "T"), c("A", "T", "C", "G"))
  testAns <- t(matrix(c(
    0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0
  ),
  nrow = 2, byrow = TRUE
  ))
  testAns <- Matrix::Matrix(testAns, sparse=  TRUE)
  #
  givenAns1 <- .one_hot_encode_sinuc(testSeqsB[[1]])
  givenAns2 <- .one_hot_encode_sinuc(testSeqsB[[2]])
  givenAns <- get_one_hot_encoded_seqs(testSeqs)
  #
  expect_is(givenAns, "dgCMatrix") ##using  Matrix
  expect_is(givenAns1, "matrix")
  expect_equal(testAns, givenAns)
})
