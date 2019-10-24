context("Auxiliary Functions II")


# test_that("get_factors_from_factor_clustering handles null hopach object", {
#     hopachObj <- NULL
#     fMat <- matrix(rep(runif(1),1000), ncol = 5)
#     expect_error(.get_factors_from_factor_clustering(hopachObj, fMat),
#                  "Hopach object is NULL")
# })


test_that("get_factors_from_factor_clustering handles all-zero
          featuresMatrix/factors", {
    hopachObj <- list(clustering = list(k = 5, sizes = rep(5,5), order = 1:25))
    fMat <- matrix(rep(0,1000), ncol = 5)
    expect_error(.get_factors_from_factor_clustering(hopachObj, fMat),
                 "WARNING: All zeroes as factors")
})


test_that("get_factors_from_factor_clustering handles NA in
          featuresMatrix/factors", {
    hopachObj <- list(clustering = list(k = 5, sizes = rep(5,5), order = 1:25))
    fMat <- matrix(rep(0,1000), ncol = 5)
    fMat[2,2] <- NA
    expect_error(.get_factors_from_factor_clustering(hopachObj, fMat),
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
