# Sys.unsetenv("R_TESTS")

context("archR main functionality")

## Make toy objects and data
inputFastaFilename <- system.file("extdata", "example_data.fa",
                                  package = "archR",
                                  mustWork = TRUE)
tssSeqs <- suppressMessages(archR::prepare_data_from_FASTA(inputFastaFilename))
tssSeqsRaw <- suppressMessages(archR::prepare_data_from_FASTA(inputFastaFilename, rawSeq = TRUE))

nSeqs <- ncol(tssSeqs)
positions <- seq(1,100)
# sinuc <- c('A', 'C', 'G', 'T')
# changedOrder <- sample.int(nSeqs, nSeqs, replace = FALSE)
# tssSeqs <- tssSeqsOriginal[ , changedOrder]
toyConfig <- archR::archRSetConfig(innerChunkSize = 100,
                            kMin = 2, kMax = 20, parallelize = FALSE,
                            modSelType = "stability", 
                            nIterationsUse = 50,
                            nCoresUse = NA)


test_that("archR works", {
    set.seed(1234)
    archRresult <- suppressMessages(archR::archR(toyConfig, seqsRaw = tssSeqsRaw, 
                                seqsMat = tssSeqs, thresholdItr = 1))
    ##
    ## We just the seqsClustLabels and a block of the basis vectors matrix
    # toyResult_seqsClustLabels <- 
    #     c("1", "1", "2", "3", "2", "1", "2", "2", "3", "2", "1", "2", "2", "3", 
    #       "3", "2", "2", "2", "2", "1", "3", "2", "1", "3", "2", "3", "3", "3", 
    #       "3", "2", "2", "1", "3", "3", "3", "2", "2", "2", "2", "3", "2", "3", 
    #       "2", "3", "2", "1", "2", "3", "2", "1", "2", "1", "2", "2", "2", "2", 
    #       "3", "2", "2", "2", "2", "3", "1", "2", "3", "3", "2", "1", "2", "2", 
    #       "1", "1", "2", "3", "1", "3", "2", "3", "3", "1", "3", "1", "1", "2", 
    #       "1", "1", "1", "1", "2", "3", "2", "1", "3", "2", "2", "3", "3", "2", 
    #       "1", "3", "3", "1", "2", "1", "1", "2", "2", "3", "2", "3", "3", "1", 
    #       "3", "1", "2", "1", "2", "2", "2", "2", "2", "2", "1", "3", "2", "2", 
    #       "3", "2", "3", "2", "1", "2", "2", "2", "2", "2", "1", "1", "3", "2", 
    #       "2", "2", "1", "2", "1", "1", "1", "1", "3", "2", "1", "2", "2", "2", 
    #       "3", "2", "3", "1", "2", "2", "2", "2", "3", "2", "1", "2", "1", "1", 
    #       "3", "3", "1", "3", "2", "3", "2", "1", "2", "3", "1", "2", "3", "1", 
    #       "2", "2", "2", "2", "3", "2", "2", "1", "2", "3", "2", "2", "2", "2", 
    #       "1", "2", "2", "2")
    # toyResult_nBasisVectors <- 3
    # toyResult_factors_10_3 <- matrix(c(0.141375993, 0.23189513, 0.31949721,
    #                             0.129934837, 0.26171279, 0.23964413,
    #                             0.162060017, 0.23019522, 0.30487274,
    #                             0.150710132, 0.18659361, 0.25073772,
    #                             0.187403520, 0.31464845, 0.26935487), 
    #                             nrow = 10, ncol = 3, byrow = TRUE)
    ##
    # expect_identical(toyResult, archRresult)
    # expect_equal(archRresult$seqsClustLabels[[1]], toyResult_seqsClustLabels)
    # expect_equal(archRresult$clustBasisVectors[[1]]$nBasisVectors, toyResult_nBasisVectors)
    # expect_equal(archRresult$clustBasisVectors[[1]]$basisVectors[1:10,], 
    #              toyResult_factors_10_3)
    expect_equal_to_reference(archRresult, "archRresult_stability_check.rds")
    ##
    ## Test cross-validation-based model selection. This needs to parallel as TRUE.
    toyConfig <- archR::archRSetConfig(innerChunkSize = 100,
                            kMin = 2, kMax = 20, parallelize = TRUE,
                            modSelType = "cv",
                            nIterationsUse = 10,
                            nCoresUse = 2,
                            flags = list(debugFlag = FALSE,
                                         verboseFlag = TRUE,
                                         plotVerboseFlag = FALSE,
                                         timeFlag = TRUE))
    set.seed(1234)
    archRresult <- archR::archR(toyConfig, seqsRaw = tssSeqsRaw,
                                seqsMat = tssSeqs, thresholdItr = 1)
    skip_on_travis()
    expect_equal_to_reference(archRresult, "archRresult_cv_check.rds")
})


test_that("Handles negative threshold iteration", {
    # toyResult <- archR(toyConfig, seqsMat = tssSeqs, thresholdItr = -1)
    expect_error(.assert_archR_thresholdIteration(-1),
                 "Expecting threshold iteration to be numeric and > 0")
    expect_error(archR(toyConfig, seqsRaw = tssSeqsRaw, seqsMat = tssSeqs, thresholdItr = -1),
                 "Expecting threshold iteration to be numeric and > 0")
})

test_that("Handles negative threshold iteration from archR main", {
    expect_error(archR(toyConfig, seqsRaw = tssSeqsRaw, seqsMat = tssSeqs, thresholdItr = -1))
    expect_error(archR(toyConfig, seqsRaw = tssSeqsRaw, seqsMat = tssSeqs, thresholdItr = -1),
                 "Expecting threshold iteration to be numeric and > 0")
})


test_that("Config handles: negative innerChunkSize", {
    expect_error(archRSetConfig(innerChunkSize = -500,
                                kMin = 2, kMax = 8, parallelize = TRUE,
                                cvFolds = 3, nIterationsUse = 50,
                                nCoresUse = 2),
                "'innerChunkSize' should be > 0")
})


test_that("Config handles: NULL flags", {
    expect_error(archRSetConfig(innerChunkSize = 500,
                                kMin = 2, kMax = 8, parallelize = TRUE,
                                cvFolds = 3, nIterationsUse = 50,
                                nCoresUse = 2,
                                flags = NULL),
                 "'flags' are NULL")
})

test_that("Config handles: improper flags", {
    expect_error(archRSetConfig(innerChunkSize = 500,
                                kMin = 2, kMax = 8, parallelize = TRUE,
                                cvFolds = 3, nIterationsUse = 50,
                                nCoresUse = 2,
                                flags = list(bugFlag = NULL,
                                             plotVerboseFlag = FALSE,
                                             verboseFlag = TRUE,
                                             timeFlag = TRUE)),
                 "Unexpected names or no elements in flags variable")
})

test_that("Config handles: improper flags", {
    expect_error(archRSetConfig(innerChunkSize = 500,
                                kMin = 2, kMax = 8, parallelize = TRUE,
                                cvFolds = 3, nIterationsUse = 50,
                                nCoresUse = 2,
                                flags = list(debugFlag = NULL,
                                             plotVerboseFlag = FALSE,
                                             verboseFlag = TRUE,
                                             timeFlag = TRUE)),
                 "Expected only LOGICAL values in flag variable,\n                        found otherwise")
})


test_that("Config handles: negative nIterations", {
    expect_error(archRSetConfig(innerChunkSize = 500,
                                kMin = 2, kMax = 8, parallelize = TRUE,
                                cvFolds = 3, nIterationsUse = -50,
                                nCoresUse = 2,
                                flags = list(debugFlag = FALSE,
                                                plotVerboseFlag = FALSE,
                                                verboseFlag = TRUE,
                                                timeFlag = TRUE)),
                 "'nIterationsUse' should be > 0")
})



test_that("Config handles: negative minSeqs", {
    expect_error(archRSetConfig(innerChunkSize = 500,
                                kMin = 2, kMax = 8, parallelize = TRUE,
                                cvFolds = 3, nIterationsUse = 50,
                                nCoresUse = 2,
                                minSeqs = -20,
                                flags = list(debugFlag = FALSE,
                                                plotVerboseFlag = FALSE,
                                                verboseFlag = TRUE,
                                                timeFlag = TRUE)),
                    "'misSeqs' should be > 0")
})


test_that("Config handles: negative alphaVal", {
    expect_error(archRSetConfig(innerChunkSize = 500,
                                kMin = 2, kMax = 8, parallelize = TRUE,
                                cvFolds = 3, nIterationsUse = 50,
                                nCoresUse = 2,
                                minSeqs = 20,
                                alphaBase = -3,
                                alphaPow = 3,
                                flags = list(debugFlag = FALSE,
                                                plotVerboseFlag = FALSE,
                                                verboseFlag = TRUE,
                                                timeFlag = TRUE)),
                    "Resulting alpha value is < 0. Check 'alphaBase' and\n                    'alphaPow'")
})


test_that("Config handles: innerChunkSize < nSeqs", {
    expect_error(.assert_archR_innerChunkSize_in_tandem(500,450),
                    "'innerChunkSize' should be <= number of input sequences")
})

test_that("Config handles: kFolds NULL", {
    expect_error(.assert_archR_kFolds_in_tandem(NULL,450),
                    "'kFolds' is NULL")
})

test_that("Config handles: kFolds is not numeric", {
    expect_error(.assert_archR_kFolds_in_tandem("5",450),
                 "'kFolds' should be numeric and > 0")
})

test_that("Config handles: kFolds <= nSeqs", {
    expect_error(.assert_archR_kFolds_in_tandem(451,450),
                 paste0("CV folds should be less than or equal to #sequences. ",
                 "Standard values: 3, 5, 10."))
})



test_that("next iteration chunks list has a 0-length entry", {
    toyList <- vector("list", 5)
    toyList <- lapply(seq_along(toyList),
                      function(x){
                          if ( x != 3) {
                              toyList[[x]] <- rep(seq(1:25),5)
                          }
                      })

    expect_error(.assert_archR_OK_for_nextIteration(toyList),
                 "Index 3 of zero length")
})

test_that("next iteration chunks list is NULL", {
    toyList <- NULL
    expect_error(.assert_archR_OK_for_nextIteration(toyList),
                 "Chunks for next iteration are NULL")
})
