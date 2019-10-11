context("One-hot encoding functions")

test_that("Report empty sequence", {
  #
  testSeqsNotOK <- c()
  expect_error(get_one_hot_encoded_seqs(testSeqsNotOK), "empty")
  expect_error(one_hot_encode_sinuc(testSeqsNotOK), "empty")
})

test_that("One-hot encoding is correctly done", {
  testSeqs <- c("GGCT", "ATCG")
  testSeqsB <- list(c("G", "G", "C", "T"), c("A", "T", "C", "G"))
  testAns <- matrix(c(
    0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0
  ),
  nrow = 2, byrow = T
  )
  #
  givenAns1 <- one_hot_encode_sinuc(testSeqsB[[1]])
  givenAns2 <- one_hot_encode_sinuc(testSeqsB[[2]])
  givenAns <- get_one_hot_encoded_seqs(testSeqs)
  #
  expect_is(givenAns, "matrix")
  expect_is(givenAns1, "matrix")
  expect_equal(testAns, givenAns)
})
