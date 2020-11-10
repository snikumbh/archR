context("Auxiliary Functions II")

test_that("setup_clustFactors return value", {
    toyMat <- matrix(rep(runif(1),1000), ncol = 5)
    returnVal <- .setup_clustFactors_for_archR_result(toyMat)

    expect_identical(names(returnVal), c("nBasisVectors", "basisVectors"))
    expect_equal(returnVal$nBasisVectors, 5)
    expect_equal(returnVal$basisVectors, toyMat)
})

test_that("decide_process_outer_chunk works fine", {

    expect_message(.decide_process_outer_chunk(25, 24, 5),
                    "Sorry, will not process this small a chunk: 24")
    expect_true(.decide_process_outer_chunk(25, 24, 5))
    expect_false(.decide_process_outer_chunk(25, 30, 5))
    expect_error(.decide_process_outer_chunk(15, 30, 4),
                    "'minSeqs' should be at least 4 times 'kFolds'")
    expect_error(.decide_process_outer_chunk(25, 0, 5),
                    "Outer chunk of size 0")
})

# test_that("handle_chunk_w_NMF works fine", {
#     toyList <- vector("list", 5)
#     toyList <- lapply(seq_along(toyList),
#                       function(x){
#                           # if ( x != 3) {
#                               toyList[[x]] <- rep(seq(1:25),5)
#                           # }
#                       })
#     toyiChunksColl <- toyList
#     fMat <- matrix(rep(runif(1),1000), ncol = 5)
#     config <- archRSetConfig(innerChunkSize = 500,
#                                 kMin = 2, kMax = 8, parallelize = FALSE,
#                                 cvFolds = 3, nIterationsUse = 50,
#                                 nCoresUse = NA, seedVal = 10208090,
#                                 flags = list(debugFlag = FALSE,
#                                              plotVerboseFlag = FALSE,
#                                              verboseFlag = TRUE,
#                                              timeFlag = TRUE))
#     config$flags$debugFlag = NULL
#     expect_error(.handle_chunk_w_NMF(3, toyiChunksColl, fMat, config),
#                  "Expected only LOGICAL values in flag variable,\n
# found otherwise")
#
# })

test_that("get_dimers is working fine", {
  alphabet = c("A", "C", "G", "T")
  expectAns <- c("AA", "CA", "GA", "TA", "AC", "CC", "GC", "TC", "AG", "CG", 
                 "GG", "TG", "AT", "CT", "GT", "TT")
  ans <- get_dimers_from_alphabet(alphabet)
  expect_equal(expectAns, ans)
  expect_error(get_dimers_from_alphabet(NULL))
})


test_that("get_hopach_cluster_medoidsIdx handles null hopach object", {
    hopachObj <- NULL
    fMat <- matrix(rep(runif(1),1000), ncol = 5)
    expect_error(.get_hopach_cluster_medoidsIdx(hopachObj),
                 "Hopach object is NULL")
    ## Make hopach object to test

})


test_that("NMF result from py script OK", {
    nmfList <- vector("list", 2)
    samplesMat <- matrix(rep(runif(1),1000), ncol = 5)
    featuresMat <- matrix(rep(runif(1),1000), nrow = 5)
    expect_error(get_features_matrix(nmfList), "0LengthEntry")
    nmfList[[1]] <- featuresMat
    expect_error(get_samples_matrix(nmfList), "0LengthEntry")
    nmfList[[2]] <- samplesMat
    expect_equal(get_samples_matrix(nmfList), samplesMat)
    expect_equal(get_features_matrix(nmfList), featuresMat)
})


test_that("getting dimers works", {
    expect_error(get_dimers_from_alphabet(NULL))
    expect_equal(get_dimers_from_alphabet(c(1,0)), c("11", "01", "10", "00"))
})




test_that("get_factors_from_factor_clustering handles null hopach object", {
    hopachObj <- NULL
    fMat <- matrix(rep(runif(1),1000), ncol = 5)
    expect_equal(.get_factors_from_factor_clustering2(hopachObj, fMat),
                 fMat)
})


test_that("get_factors_from_factor_clustering handles all-zero
          featuresMatrix/factors", {
    hopachObj <- list(clustering = list(k = 5, sizes = rep(5,5), order = 1:25))
    fMat <- matrix(rep(0,1000), ncol = 5)
    expect_error(.get_factors_from_factor_clustering2(hopachObj, fMat),
                 "WARNING: All zeroes as factors")
})


test_that("get_factors_from_factor_clustering handles NA in
          featuresMatrix/factors", {
    hopachObj <- list(clustering = list(k = 5, sizes = rep(5,5), order = 1:25))
    fMat <- matrix(rep(0,1000), ncol = 5)
    fMat[2,2] <- NA
    expect_error(.get_factors_from_factor_clustering2(hopachObj, fMat),
                 "Factors have NA")
})


# test_that("handle_chunk_w_NMF handles invalid nCores", {
#     innerChunkIdx <- -1
#     tempList <- vector("list", 5)
#     innerChunksColl <- lapply(seq_along(tempList),
#                               function(x){
#                                   tempList[[x]] <- round(10*runif(5))
#                               })
#     this_mat <- matrix(rep(0,1000), ncol = 5)
#     config <- archRSetConfig(parallelize = TRUE, nCoresUse = -32)
#     expect_error(.handle_chunk_w_NMF(innerChunkIdx, innerChunksColl,
#                                      this_mat, config),
#                  "Specified more than available cores. Stopping ")
# })
