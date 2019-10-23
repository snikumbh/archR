context("Auxiliary Functions I")


# test_that("collate_clusters handles null hopach object", {
#     hopachObj <- NULL
#     tempList <- vector("list", 5)
#     globClustAssignments <- lapply(tempList,
#                                     function(x){
#                                         x <- round(10*runif(5))
#                                     })
#     expect_error(.collate_clusters(hopachObj, globClustAssignments),
#                  "Hopach object is NULL")
# })

test_that("collate_clusters handles empty globClustAssignments", {
    hopachObj <- list(clustering = list(k = 5, sizes = rep(5,5), order = 1:25))
    tempList <- vector("list", 5)
    globClustAssignments <- lapply(seq_along(tempList),
                                   function(x){
                                       if ( x != 3) {
                                            tempList[[x]] <- round(10*runif(5))
                                       }
                                   })
    expect_error(.collate_clusters(hopachObj, globClustAssignments),
                 "Cluster assignments variable has a 0-length entry")
})


test_that("update_cluster_labels handles empty collatedClustAssignments", {
    seqsClustLabels <- rep("0-1-2-3", 25)
    tempList <- vector("list", 5)
    globClustAssignments <- lapply(seq_along(tempList),
                                   function(x){
                                       if ( x != 3) {
                                           tempList[[x]] <- round(10*runif(5))
                                       }
                                   })
    expect_error(.update_cluster_labels(seqsClustLabels, globClustAssignments),
                 "Cluster assignments variable has a 0-length entry")
})


# test_that("update_cluster_labels handles inconsistent #sequences", {
#     seqsClustLabels <- rep("0-1-2-3", 20)
#     tempList <- vector("list", 5)
#     globClustAssignments <- lapply(seq_along(tempList),
#                                    function(x){
#                                     tempList[[x]] <- round(10*runif(5))
#                                    })
#     expect_error(.update_cluster_labels(seqsClustLabels, globClustAssignments),
#     "Number of sequences in seqsClustLabels and clustAssignments not equal")
# })


test_that("prepare_chunks handles negative chunkSize", {
    seqsClustLabels <- rep("0-1-2-3", 20)
    tempList <- vector("list", 5)
    globClustAssignments <- lapply(seq_along(tempList),
                                   function(x){
                                       tempList[[x]] <- round(10*runif(5))
                                   })
    expect_error(.prepare_chunks(seqsClustLabels, -25),
                 "'innerChunkSize' should be > 0")
})





test_that("handle_clustering_of_factors handles all-zero featuresMatrix/factors", {

    fMat <- matrix(rep(0,1000), ncol = 5)
    expect_error(.handle_clustering_of_factors(fMat),
                   "WARNING: All zeroes as factors")
})


test_that("handle_clustering_of_factors handles NA in featuresMatrix/factors", {

    fMat <- matrix(rep(0,1000), ncol = 5)
    fMat[2,2] <- NA
    expect_error(.handle_clustering_of_factors(fMat),
                   "Factors have NA")
})

test_that("handle_clustering_of_factors handles improper flags", {

    fMat <- matrix(rep(0,1000), ncol = 5)
    fMat[2,2] <- NA
    expect_error(.handle_clustering_of_factors(fMat, flags = list()),
                 "")
})




#
# test_that("Q^2 computation: handling original non-matrix", {
#     matA <- c(rnorm(20)) # matrix(rnorm(20), nrow = 4)#rnorm(20), nrow = 4)
#     reMatA <- matrix(rnorm(20), nrow = 4)
#     expect_error(.compute_q2(matA, reMatA), "not of type matrix")
# })
#
# test_that("Q^2 computation: handling reconstructed non-matrix", {
#     matA <- matrix(rnorm(20), nrow = 4)
#     reMatA <- c(rnorm(20)) # matrix()
#     expect_error(.compute_q2(matA, reMatA), "not of type matrix")
# })

# test_that("NMF solved correctly", {
#
#
#
# })
