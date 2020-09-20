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
    # set.seed(1234)
    # toyClustLabels <- sample(x = rep(as.character(1:25),3), size = 75,
    #                          replace = FALSE)
    toyClustLabels <- c("9", "22", "20", "23", "12", "22", "1", "16", "20",
                        "9", "21", "10", "18", "8", "13", "1", "17", "18",
                        "11", "14", "11", "9", "25", "3", "12", "16", "1",
                        "19", "15", "2", "21", "10", "6", "22", "8", "6",
                        "25", "10", "12", "5", "20", "17", "7", "14", "8",
                        "21", "7", "18", "7", "4", "2", "14", "4", "23", "4",
                        "24", "13", "3", "15", "5", "13", "19", "5", "11", "23",
                        "24", "24", "15", "16", "17", "2", "6", "19", "25", "3")
    # collection <- sample.int(25)
    collection <- c("3", "15", "24", "14", "19", "13", "1", "5", "12", "9",
                    "11", "8", "4", "17", "20", "16", "25", "10", "2", "7",
                    "6", "18", "21", "22", "23")
    tempList <- vector("list", 5)
    globClustAssignmentsErr <- lapply(seq_along(tempList),
                                   function(x){
                                       if ( x != 3) {
                                           tempList[[x]] <- round(10*runif(5))
                                       }
                                   })
    # set.seed(1234)
    # collection <- sample.int(25)
    idx <- rep(1:5,5)
    globClustAssignments <- lapply(seq_along(tempList),
                                   function(x){
                                       tempList[[x]] <- collection[idx == x]
                                   })
    ## this tests an assertion
    expect_error(.update_cluster_labels(toyClustLabels,
                                        globClustAssignmentsErr),
                 "Cluster assignments variable has a 0-length entry")
    ####
    #### Test flag set to NULL instead of expected logical
    toyFlags <- list(debugFlag = NULL,
                        plotVerboseFlag = FALSE,
                        verboseFlag = FALSE,
                        timeFlag = TRUE)

    ## this tests the flags' assertion
    expect_error(.update_cluster_labels(toyClustLabels, globClustAssignments,
                                        toyFlags),
                    "Expected only LOGICAL values in flag variable,
                        found otherwise")
    #### Setting flags properly
    toyFlags$debugFlag <- FALSE
    toyFlags$verboseFlag <- FALSE

    samarth <- c("5", "4", "5", "5", "4", "4", "2", "1", "5", "5", "3", "3",
                 "2", "2", "1", "2", "4", "2", "1", "4", "1", "5", "2", "1",
                 "4", "1", "2", "5", "2", "4", "3", "3", "1", "4", "2", "1",
                 "2", "3", "4", "3", "5", "4", "5", "4", "2", "3", "5", "2",
                 "5", "3", "4", "4", "3", "5", "3", "3", "1", "1", "2", "3",
                 "1", "5", "3", "1", "5", "3", "3", "2", "1", "4", "4", "1",
                 "5", "2", "1")
    ## this tests the updated cluster labels
    ## -- the updated labels should have 5 unique cluster labels since
    ## globClustAssignments has 5 clusters
    ## --
    print("\n")
    print(samarth)
    print("\n")
    samarth_ans <- .update_cluster_labels(toyClustLabels, globClustAssignments,
                                          toyFlags)
    # cat(paste(shQuote(samarth_ans, type="cmd"), collapse=", "))

    expect_equal(samarth_ans, samarth)

    expect_identical(sort(unique(samarth_ans)), sort(unique(samarth)))
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



test_that("get seq clusters as list works fine", {
    set.seed(1234)
    toyClustLabels <- sample(x = rep(as.character(1:14),3), size = 42, replace = FALSE)
    seqsClustLabels <- toyClustLabels
    expAns <- 14
    expect_equal(length(get_seqs_clusters_in_a_list(seqsClustLabels)), 14)
})




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





test_that("clustering of factors handles all-zero featuresMatrix", {

    fMat <- matrix(rep(0,1000), ncol = 5)
    expect_error(.handle_clustering_of_factors(fMat),
                   "WARNING: All zeroes as factors")
})


test_that("handle_clustering_of_factors handles NA in featuresMatrix", {

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
