context("Representing sequences as an image")

## TODO: Add a visual test using the new func viz_seqs_as_acgt_mat_from_seqs

test_that("Handling image representation of sequences", {
    testSeqs_filename <-
        system.file("extdata",
                    "example_data.fa",
                    package = "archR",
                    mustWork = TRUE)
    testSeqsOK <- prepare_data_from_FASTA(testSeqs_filename, rawSeq = TRUE)
    positionLabelsTest <- seq(1, 100)
    vdiffr::expect_doppelganger(
        "Sample image of sequences",
        viz_seqs_as_acgt_mat_from_seqs(as.character(testSeqsOK),
                                position_labels = positionLabelsTest,
                                xt_freq = 25, yt_freq = 200,
                                savefilename = NULL)
    )
})

# test_that("Handling 0/negative xt_freq", {
#     testSeqs_filename <-
#         system.file("extdata",
#             "example_data.fa",
#             package = "archR",
#             mustWork = TRUE)
#     testSeqsOK <- prepare_data_from_FASTA(testSeqs_filename, rawSeq = TRUE)
#     positionLabelsTest <- seq(1, 100)
#     expect_warning(viz_seqs_as_acgt_mat_from_seqs(as.character(testSeqsOK),
#                         position_labels = positionLabelsTest,
#                         xt_freq = -5, yt_freq = 100,
#                         savefilename = NULL),
#         regexp = "Expected positive integer (< length of sequences) for"
#     )
# })
# 
# test_that("Handling 0/negative yt_freq", {
#     testSeqs_filename <-
#         system.file("extdata",
#             "example_data.fa",
#             package = "archR",
#             mustWork = TRUE)
#     testSeqsOK <- prepare_data_from_FASTA(testSeqs_filename, rawSeq = TRUE)
#     positionLabelsTest <- seq(1, 100)
#     expect_warning(viz_seqs_as_acgt_mat_from_seqs(as.character(testSeqsOK),
#                         position_labels = positionLabelsTest,
#                         xt_freq = 25, yt_freq = 0,
#                         savefilename = NULL),
#         regexp = "Expected positive integer (< number of sequences) for"
#         )
# })


