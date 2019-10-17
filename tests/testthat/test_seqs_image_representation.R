context("Representing sequences as an image")

# test_that("Handling image representation of sequences", {
#     testSeqs_filename <-
#         system.file("extdata",
#                     "example_data.fa",
#                     package = "archR",
#                     mustWork = TRUE)
#     testSeqsOK <- prepare_data_from_FASTA(testSeqs_filename, rawSeq = TRUE)
#     positionLabelsTest <- setdiff(seq(-100, 100), c(0))
#     vdiffr::expect_doppelganger(
#         "Sample image of sequences",
#         represent_matrix_of_acgt_with_ggplot2(as.character(testSeqsOK),
#                                  position_labels = positionLabelsTest,
#                                  annClusters = NULL,
#                                  saveFilename = NULL)
#     )
# })
