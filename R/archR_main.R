#' @title
#' This function processes the given data set.
#'
#' @description Set archR configuration and call this function to process a
#' data set using archR.
#'
#' @param config archR configuration object
#' @param tss.seqs sequences
#' @param thresholdItr number of iterations to perform
#'
#' @return Object of class archR (or at the moment, a list of lists)
#' @export
archR <- function(config, tss.seqs, thresholdItr = 2) {
    ### Assume: config should have a valid, default set of params
    ### Make checks in the wrapper functions
    ### To continue archR from an earlier run (further levels
    ### downstream)
    ### 1. Initializations of seqClustLabels, Factors etc should be
    ### appropriately handled.
    ### Initialize sequence-cluster-labels and outerChunks
    seqsClustLabels <- rep("0", ncol(tss.seqs))
    clustFactors <- vector("list", thresholdItr)
    # intClustFactors <- NULL
    outerChunksColl <- vector("list", 1)
    ### Set outerChunks for first iteration
    outerChunksColl[[1]] <- seq(ncol(tss.seqs))
    test_itr <- 0
    ##--------------------------------------------------------------------------
    while (test_itr < thresholdItr) {
        message("=== Processing Level ", test_itr, ", ",
                length(outerChunksColl), " chunk(s) ===")
        nxtOuterChunksColl <- vector("list")
        intClustFactors <- NULL
        for (outerChunkIdx in 1:length(outerChunksColl)) {
            ###
            ###
            outerChunk <- outerChunksColl[[outerChunkIdx]]
            doNotProcess <- decide_process_outer_chunk(config$minSeqs,
                                                       length(outerChunk),
                                                       config$kFolds)
            ###
            if (doNotProcess) {
                collatedClustAssignments <- list(outerChunk)
                # print("=== dim of clustFactors prev itr ===")
                # print(dim(clustFactors[[test_itr]]$Factors))
                if (!is.null(intClustFactors)) {
                  intClustFactors <-
                      cbind(intClustFactors,
                            as.matrix(
                                clustFactors[[test_itr]]$Factors[,outerChunkIdx]))
                } else {
                  intClustFactors <-
                      as.matrix(
                          clustFactors[[test_itr]]$Factors[, outerChunkIdx]
                          )
                }

            } else {
                innerChunksColl <- prepare_chunks(outerChunk,
                                                  config$innerChunkSize,
                                                  length(outerChunk))
                ###
                globFactors <- vector("list", length(innerChunksColl))
                globClustAssignments <- vector("list", length(innerChunksColl))
                ###
                for (innerChunkIdx in 1:length(innerChunksColl)) {
                  # print("=== inner chunk size: ===")
                  # print(length(innerChunksColl[[innerChunkIdx]]))
                  this_tss.seqs <- tss.seqs[, innerChunksColl[[innerChunkIdx]]]
                  thisNMFResult <- handle_chunk_w_NMF(innerChunkIdx,
                                                      innerChunksColl,
                                                      this_tss.seqs,
                                                      globFactors,
                                                      globClustAssignments,
                                                      config)
                  globFactors <- thisNMFResult$globFactors
                  globClustAssignments <- thisNMFResult$globClustAssignments
                }
                ### for loop over innerChunksColl ENDS
                ### need globFactors, globClustAssignments
                globClustAssignments <- unlist(globClustAssignments,
                                               recursive = FALSE)
                ###
                ### CBind the factors from all inner chunks into one matrix
                globFactorsMat <- do.call(cbind, globFactors)
                ###
                hopachDecision <- decide_hopach(globFactorsMat,
                                                distMethod = "cosangle",
                                                withinMeasure = "mean")
                ###
                globFactorsClustering <- NULL
                if (hopachDecision) {
                  ### We need to cluster the factors
                  globFactorsClustering <-
                      handle_clustering_of_factors(globFactorsMat,
                                                   distMethod = "cosangle",
                                                   flags = config$flags)
                }
                ### Manage factors
                if (!is.null(intClustFactors)) {
                  intClustFactors <-
                      cbind(intClustFactors,
                            get_factors_from_factor_clustering(globFactorsClustering,
                                globFactorsMat))
                } else {
                  intClustFactors <-
                      get_factors_from_factor_clustering(globFactorsClustering,
                                                         globFactorsMat)
                }
                ### Manage collated cluster assignments
                collatedClustAssignments <- collate_clusters(globFactorsClustering,
                  globClustAssignments)
                ###
            }  ### IfElse doNotProcess outer chunk ENDS
            ###
            seqsClustLabels <-
                update_cluster_labels(seqsClustLabels,
                                      collatedClustAssignments =
                                          collatedClustAssignments,
                                      flags = config$flags)

            ### Collect (append) clusters at current level
            nxtOuterChunksColl <- append(nxtOuterChunksColl,
                                         collatedClustAssignments)
            message(paste0("Outer chunk ", outerChunkIdx, ", current total basis vectors: ",
                           ncol(intClustFactors)))

        }  ### for loop over outerChunksCollection ENDS
        #
        clustFactors[[test_itr + 1]] <-
            list(nBasisVectors = ncol(intClustFactors),
                 basisVectors = intClustFactors)
        outerChunksColl <- nxtOuterChunksColl
        test_itr <- test_itr + 1
    }  # algorithm while loop ENDS

    archRresult <- list(seqsClustLabels = seqsClustLabels, clustBasisVectors = clustFactors,
        config = config, call = match.call())
    message("=== archR exiting, returning result ===")
    return(archRresult)
}  #archR function ENDS
