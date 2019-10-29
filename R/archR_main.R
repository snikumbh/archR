#' @title
#' This function processes the given data set.
#'
#' @description Call this function to process a data set using archR.
#'
#' @param config archR configuration object
#' @param seqsMat A matrix. A matrix of one-hot encoded sequences with sequences
#' along columns
#' @param thresholdItr numeric. Specify the number of iterations to perform.
#' Default is 2
#'
#' @return Object of class archR (or at the moment, a list of lists)
#' @export
archR <- function(config, seqsMat, thresholdItr = 2) {
    ## Make checks for params in configuration
    .assert_archR_config(config, ncol(seqsMat))
    ## ** To continue archR from an earlier run (further levels downstream)
    ## 1. Initializations of seqClustLabels, Factors etc should be
    ## appropriately handled.**
    ## Initialize sequence-cluster-labels and outerChunks
    seqsClustLabels <- rep("0", ncol(seqsMat))
    clustFactors <- vector("list", thresholdItr)
    outerChunksColl <- vector("list", 1)
    ## Set outerChunks for first iteration
    outerChunksColl[[1]] <- seq(ncol(seqsMat))
    test_itr <- 0
    ##--------------------------------------------------------------------------
    while (test_itr < thresholdItr) {
        message("=== Processing Level ", test_itr, ", ",
                length(outerChunksColl), " chunk(s) ===")
        nxtOuterChunksColl <- vector("list")
        intClustFactors <- NULL
        for (outerChunkIdx in seq_along(outerChunksColl)) {
            ##
            outerChunk <- outerChunksColl[[outerChunkIdx]]
            ## Make a decision to process based on size of chunk
            message("Deciding whether to process chunk?")
            doNotProcess <- .decide_process_outer_chunk(config$minSeqs,
                                                        length(outerChunk),
                                                        config$kFolds)
            ##
            if (doNotProcess) {
                ## TO-DO: Could write function to manipulate
                ## collatedClustAssignments
                ## NOTE: We enter here only from second iteration onwards,
                ## otherwise statements here can fail. Because, clustFactors is
                ## assumed to be populated already -- in the first iteration,
                ## test_itr is 0, and indexing using test_itr would fail.
                collatedClustAssignments <- list(outerChunk)
                if (!is.null(intClustFactors)) {
                    intClustFactors <- cbind(intClustFactors,
                        as.matrix(
                        clustFactors[[test_itr]]$basisVectors[, outerChunkIdx])
                        )
                } else {
                    intClustFactors <- as.matrix(
                        clustFactors[[test_itr]]$basisVectors[, outerChunkIdx]
                        )
                }
            } else {
                message("Yes")
                innerChunksColl <- .prepare_chunks(outerChunk,
                                                    config$innerChunkSize)
                ## Maintain these in a list, for collation later when
                ## all innerChunks in innerChunksColl have been processed
                globFactors <- vector("list", length(innerChunksColl))
                globClustAssignments <- vector("list", length(innerChunksColl))
                ##
                for (innerChunkIdx in seq_along(innerChunksColl)) {
                    ## Setting up sequences for the current chunk
                    this_seqsMat <-
                        seqsMat[, innerChunksColl[[innerChunkIdx]]]
                    thisNMFResult <- .handle_chunk_w_NMF(innerChunkIdx,
                                                        innerChunksColl,
                                                        this_seqsMat,
                                                        config)
                    .assert_archR_NMFresult(thisNMFResult)
                    globFactors[[innerChunkIdx]] <- thisNMFResult$forGlobFactors
                    globClustAssignments[[innerChunkIdx]] <-
                        thisNMFResult$forGlobClustAssignments
                }
                ## for loop over innerChunksColl ENDS
                ## We need globFactors, globClustAssignments
                ## Single unlist of globClustAssignments brings together
                ## clusters from different innerChunks into one collection
                globClustAssignments <- unlist(globClustAssignments,
                                                recursive = FALSE)
                ## CBind the factors from all inner chunks into one matrix
                globFactorsMat <- do.call(cbind, globFactors)
                ## These factors collected from all innerChunks may need
                ## clustering
                hopachDecision <- .decide_hopach(globFactorsMat,
                                                distMethod = "cosangle",
                                                withinMeasure = "mean")
                ##
                globFactorsClustering <- NULL
                if (hopachDecision) {
                    ## Cluster the factors using hopach
                    globFactorsClustering <-
                        .handle_clustering_of_factors(globFactorsMat,
                                                    distMethod = "cosangle",
                                                    flags = config$flags)
                }
                ## Manage factors
                if (!is.null(intClustFactors)) {
                    intClustFactors <-
                        cbind(intClustFactors,
                            .get_factors_from_factor_clustering(
                                globFactorsClustering,
                                globFactorsMat))
                } else {
                    intClustFactors <-
                    .get_factors_from_factor_clustering(globFactorsClustering,
                                                            globFactorsMat)
                }
                ## Manage collated cluster assignments
                collatedClustAssignments <-
                    .collate_clusters(globFactorsClustering,
                                        globClustAssignments)
                ##
            }  ## IfElse doNotProcess outer chunk ENDS
            ##
            .assert_archR_globClustAssignments(collatedClustAssignments)
            seqsClustLabels <- .update_cluster_labels(seqsClustLabels,
                                                collatedClustAssignments,
                                                flags = config$flags)
            ##
            ## Collect (append) clusters at current level
            nxtOuterChunksColl <- append(nxtOuterChunksColl,
                                        collatedClustAssignments)
            message("Outer chunk ", outerChunkIdx,
            " done, \ncurrent total basis vectors: ", ncol(intClustFactors),
            "\ncurrent total chunks for next iteration: ",
            length(nxtOuterChunksColl), "\n")
        }  ## for loop over outerChunksCollection ENDS
        ##
        clustFactors[[test_itr + 1]] <-
            .setup_clustFactors_for_archR_result(intClustFactors)
        ##
        .assert_archR_OK_for_nextIteration(nxtOuterChunksColl)
        outerChunksColl <- nxtOuterChunksColl
        test_itr <- test_itr + 1
    }  ## algorithm while loop ENDS
    ##
    temp_archRresult <- list(seqsClustLabels = seqsClustLabels,
                        clustBasisVectors = clustFactors,
                        config = config, call = match.call())
    ## reordering returns the result object with an additional field
    message("Reordering archR clusters")
    archRresult <- reorder_archRresult(temp_archRresult)
    ##
    message("=== archR exiting, returning result ===")
    return(archRresult)
}  ## archR function ENDS
## =============================================================================
