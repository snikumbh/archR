# Sys.unsetenv("R_TESTS")

context("archR main functionality")

# library(archR)
# reticulate::source_python(system.file("python/perform_nmf.py",
#                                       package = "archR",
#                                       mustWork = TRUE)
# )
# seed_val <- 11992288
# set.seed(seed_val)
## Make toy objects and data
inputFastaFilename <- system.file("extdata", "example_data.fa",
                                  package = "archR",
                                  mustWork = TRUE)
tssSeqsOriginal <- archR::prepare_data_from_FASTA(inputFastaFilename)

nSeqs <- ncol(tssSeqsOriginal)
positions <- seq(1,200)
sinuc <- c('A', 'C', 'G', 'T')
changedOrder <- sample.int(nSeqs, nSeqs, replace = FALSE)
tssSeqs <- tssSeqsOriginal[ , changedOrder]
toyConfig <- archRSetConfig(innerChunkSize = 500,
                            kMin = 2, kMax = 5, parallelize = FALSE,
                            cvFolds = 3, nIterationsUse = 5,
                            nCoresUse = NA, seedVal = 10208090)


# test_that("archR works", {
#     toyResult <- archR::archR(toyConfig, seqsMat = tssSeqs, thresholdItr = 2)
#     expect_identical(toyResult, archRresult)
# })


test_that("Handles negative threshold iteration", {
    # toyResult <- archR(toyConfig, seqsMat = tssSeqs, thresholdItr = -1)
    expect_error(.assert_archR_thresholdIteration(-1),
                 "Expecting threshold iteration to be numeric and > 0")
    expect_error(archR(toyConfig, seqsMat = tssSeqs, thresholdItr = -1),
                 "Expecting threshold iteration to be numeric and > 0")
})

# test_that("Handles negative threshold iteration from archR main", {
#     expect_error(archR(toyConfig, seqsMat = tssSeqs, thresholdItr = -1))
#     # expect_error(archR(toyConfig, seqsMat = tssSeqs, thresholdItr = -1),
#     #              "Expecting threshold iteration to be numeric and > 0")
# })


test_that("Config handles: negative innerChunkSize", {
    expect_error(archRSetConfig(innerChunkSize = -500,
                                kMin = 2, kMax = 8, parallelize = TRUE,
                                cvFolds = 3, nIterationsUse = 50,
                                nCoresUse = 2, seedVal = 10208090),
                "'innerChunkSize' should be > 0")
})


test_that("Config handles: NULL flags", {
    expect_error(archRSetConfig(innerChunkSize = 500,
                                kMin = 2, kMax = 8, parallelize = TRUE,
                                cvFolds = 3, nIterationsUse = 50,
                                nCoresUse = 2, seedVal = 10208090,
                                flags = NULL),
                 "'flags' are NULL")
})

test_that("Config handles: improper flags", {
    expect_error(archRSetConfig(innerChunkSize = 500,
                                kMin = 2, kMax = 8, parallelize = TRUE,
                                cvFolds = 3, nIterationsUse = 50,
                                nCoresUse = 2, seedVal = 10208090,
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
                                nCoresUse = 2, seedVal = 10208090,
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
                                nCoresUse = 2, seedVal = 10208090,
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
                                nCoresUse = 2, seedVal = 10208090,
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
                                nCoresUse = 2, seedVal = 10208090,
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
