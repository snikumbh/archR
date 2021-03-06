context("Representing sequences as an image")

## TODO: Add a visual test using the new func viz_seqs_as_acgt_mat_from_seqs

test_that("Handling image representation of sequences", {
    testSeqs_filename <- system.file("extdata",
                            "example_data.fa",
                            package = "archR",
                            mustWork = TRUE)
    testSeqsOK <- suppressMessages(archR::prepare_data_from_FASTA(
                            testSeqs_filename, raw_seq = TRUE))
    positionLabelsTest <- seq(1, 100)
    vdiffr::expect_doppelganger(
            "Sample image of sequences",
            suppressMessages(
                archR::viz_seqs_acgt_mat_from_seqs(as.character(testSeqsOK),
                                    pos_lab = positionLabelsTest,
                                    xt_freq = 25, yt_freq = 200))
    )
})
