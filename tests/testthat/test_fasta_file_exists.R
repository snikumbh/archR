context("Handling file (non-)existance")

test_that("File non-existance handling", {

  testSeqs_filenameNotOK <- system.file("extdata", "example_data_.fa")
  # This returns empty string when file not found

  expect_error(prepare_data_from_FASTA(testSeqs_filenameNotOK), "not found")

})
