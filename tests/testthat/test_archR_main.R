# Sys.unsetenv("R_TESTS")


context("archR main functionality")
library(archR)




test_that("archR (stability) works when timeFlag is FALSE", {
    ## Make toy objects and data
    inputFastaFilename <- system.file("extdata", "example_data.fa",
                                      package = "archR",
                                      mustWork = TRUE)
    
    
    
    tssSeqs_sinuc <- suppressMessages(archR::prepare_data_from_FASTA(inputFastaFilename))
    tssSeqsRaw <- suppressMessages(archR::prepare_data_from_FASTA(inputFastaFilename, rawSeq = TRUE))
    
    nSeqs <- ncol(tssSeqs_sinuc)
    positions <- seq(1,100)
    
    useFlags <- list(debugFlag = FALSE,
                     verboseFlag = TRUE,
                     plotVerboseFlag = FALSE,
                     timeFlag = FALSE)
    toyConfig <- archR::archRSetConfig(innerChunkSize = 100,
                                       kMin = 2, kMax = 20, parallelize = FALSE,
                                       modSelType = "stability", 
                                       nIterationsUse = 50,
                                       nCoresUse = NA,
                                       flags = useFlags)
    set.seed(1234)
    archRresult <- suppressMessages(archR::archR(toyConfig, seqsRaw = tssSeqsRaw, 
                                seqsMat = tssSeqs_sinuc, thresholdItr = 1))
    ##
    expect_equal_to_reference(archRresult, "archRresult_stability_check_timeFalse.rds")
    ##
})


test_that("archR (cv) works when timeFlag is FALSE", {
    ## Make toy objects and data
    inputFastaFilename <- system.file("extdata", "example_data.fa",
                                      package = "archR",
                                      mustWork = TRUE)
    
    
    
    tssSeqs_sinuc <- suppressMessages(archR::prepare_data_from_FASTA(inputFastaFilename))
    tssSeqsRaw <- suppressMessages(archR::prepare_data_from_FASTA(inputFastaFilename, rawSeq = TRUE))
    
    nSeqs <- ncol(tssSeqs_sinuc)
    positions <- seq(1,100)
    
    useFlags <- list(debugFlag = FALSE,
                     verboseFlag = TRUE,
                     plotVerboseFlag = FALSE,
                     timeFlag = FALSE)
    toyConfig <- archR::archRSetConfig(innerChunkSize = 100,
                                       kMin = 2, kMax = 20, parallelize = TRUE,
                                       modSelType = "cv",
                                       nIterationsUse = 10,
                                       nCoresUse = 2,
                                       flags = useFlags)
    ## Test cross-validation-based model selection. This needs to parallel as TRUE.
    # skip_on_travis()
    set.seed(1234)
    archRresult <- suppressMessages(archR::archR(toyConfig, seqsRaw = tssSeqsRaw,
                                        seqsMat = tssSeqs_sinuc, thresholdItr = 1))
    expect_equal_to_reference(archRresult, "archRresult_cv_check_timeFalse.rds")
})



test_that("archR (stability) works when debug & timeFlag is FALSE", {
    ## Make toy objects and data
    inputFastaFilename <- system.file("extdata", "example_data.fa",
                                      package = "archR",
                                      mustWork = TRUE)
    
    
    
    tssSeqs_sinuc <- suppressMessages(archR::prepare_data_from_FASTA(inputFastaFilename))
    tssSeqsRaw <- suppressMessages(archR::prepare_data_from_FASTA(inputFastaFilename, rawSeq = TRUE))
    
    nSeqs <- ncol(tssSeqs_sinuc)
    positions <- seq(1,100)
    ## keeping debugFlag as TRUE doesn't alter the archR result object which we 
    ## check, but it helps hit many more debug message lines in the code
    useFlags <- list(debugFlag = TRUE,
                     verboseFlag = TRUE,
                     plotVerboseFlag = FALSE,
                     timeFlag = FALSE)
    toyConfig <- archR::archRSetConfig(innerChunkSize = 100,
                                       kMin = 2, kMax = 20, parallelize = FALSE,
                                       modSelType = "stability", 
                                       nIterationsUse = 50,
                                       nCoresUse = NA,
                                       flags = useFlags)
    set.seed(1234)
    archRresult <- suppressMessages(archR::archR(toyConfig, seqsRaw = tssSeqsRaw, 
                                    seqsMat = tssSeqs_sinuc, thresholdItr = 1))
    ##
    expect_equal_to_reference(archRresult, "archRresult_stability_check_debugTrue.rds")
    ##
})


test_that("Handles negative threshold iteration", {
    # toyResult <- archR(toyConfig, seqsMat = tssSeqs, thresholdItr = -1)
    expect_error(.assert_archR_thresholdIteration(-1),
                 "Expecting threshold iteration to be numeric and > 0")
    expect_error(archR(toyConfig, seqsRaw = tssSeqsRaw, seqsMat = tssSeqs_sinuc, thresholdItr = -1),
                 "Expecting threshold iteration to be numeric and > 0")
})

test_that("Handles negative threshold iteration from archR main", {
    expect_error(archR(toyConfig, seqsRaw = tssSeqsRaw, seqsMat = tssSeqs_sinuc, thresholdItr = -1))
    expect_error(archR(toyConfig, seqsRaw = tssSeqsRaw, seqsMat = tssSeqs_sinuc, thresholdItr = -1),
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
