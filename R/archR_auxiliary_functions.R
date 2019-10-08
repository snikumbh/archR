collect_cluster_labels <- function(given_seqsClustLabels, choose_levels = 1) {
    ### Check if all_ok, all elements should have same length
    splitChar <- "-"
    elements_length <- unique(unlist(lapply(strsplit(given_seqsClustLabels,
                                                     split = splitChar),
        length)))

    assertthat::are_equal(length(elements_length), 1)

    if (choose_levels > elements_length) {
        stop(paste0("choose_levels(", choose_levels, ") greater than levels
                    present in cluster_labels(", elements_length, ")."))
    } else {
        selectedLabels <- unlist(lapply(strsplit(given_seqsClustLabels,
                                                 split = splitChar),
            function(x) {
                paste0(x[1:choose_levels], collapse = splitChar)
            }))
    }
    return(selectedLabels)
}

## =============================================================================

map_clusters_to_factors <- function(samplesMatrix, clustOrderIdx, iChunksColl,
                                    iChunkIdx, flags) {
    # clustOrderIdx IS clustering_sol_kmeans$reordering_idx the reorderingIdx is a
    # nested list, i.e. a list holding list of seqs_ids belonging to one cluster
    # together. Hence, length(this variable) gives #clusters samplesMatrix IS
    # samplesMatrix
    if (length(clustOrderIdx) == 1) {

        rightClusterOrders <- vector("list", length(clustOrderIdx))

        rightClusterOrders[[1]] <- iChunksColl[[iChunkIdx]][clustOrderIdx[[1]]]

        return(rightClusterOrders)

    } else if (length(clustOrderIdx) > 1) {
        rightClusterOrders <- vector("list", length(clustOrderIdx))
        for (cluster_idx in 1:length(clustOrderIdx)) {
            relevant_factor <-
                which.max(rowMeans(
                as.matrix(samplesMatrix[, clustOrderIdx[[cluster_idx]]])
                ))
            if (flags$debugFlag) {
                print(relevant_factor)
            }
            if (length(relevant_factor) > 1) {
                if (flags$debugFlag) {
                  print("Factor-mapping NOT OK")
                }
            } else {
                if (flags$debugFlag) {
                  print("Factor-mapping OK")
                }
                rightClusterOrders[[relevant_factor]] <-
                    iChunksColl[[iChunkIdx]][clustOrderIdx[[cluster_idx]]]
            }
        }

        return(rightClusterOrders)
    }
}
## =============================================================================

collate_clusters <- function(hopachObj, globClustAssignments) {
    ### - Based on which factors are being combined (hopachObj), look at the
    ### globClustAssignments variable to combine the respective sequences
    ### together.

    if (is.null(hopachObj)) {
        # When HOPACH not performed, the global_cluster_assingments variable
        # holds a list of lists, so we should flatten it with unlist()
        collatedClustAssignments <- globClustAssignments
        # collatedClustAssignments <- unlist(globClustAssignments,
        # recursive = FALSE)

    } else {
        ###
        nClusters <- hopachObj$clustering$k
        clustSizes <- hopachObj$clustering$sizes
        elementsOrder <- hopachObj$clustering$order
        ###
        collatedClustAssignments <- vector("list", nClusters)
        ###
        for (coll_cluster_idx in 1:length(collatedClustAssignments)) {
            to_combine_this_itr <- clustSizes[coll_cluster_idx]
            leave_first <- 0
            if (coll_cluster_idx > 1) {
                leave_first <- sum(clustSizes[1:(coll_cluster_idx - 1)])
            }
            pick_this_itr <- elementsOrder[
                (leave_first + 1):(leave_first + to_combine_this_itr)]
            #
            if (all(pick_this_itr <= length(globClustAssignments))) {
                temp <- unlist(lapply(pick_this_itr, function(x) {
                  globClustAssignments[[x]]
                }))
                collatedClustAssignments[[coll_cluster_idx]] <- temp
            } else {
                pick_this_itr <- pick_this_itr[]
                these_not_ok <-
                    which(pick_this_itr > length(globClustAssignments))
                resolve_these <- these_not_ok
                collatedClustAssignments[[coll_cluster_idx]] <-
                    unlist(lapply(pick_this_itr[-these_not_ok],
                  function(x) {
                    globClustAssignments[[x]]
                  }))
            }
        }  # collation for loop ends
    }
    return(collatedClustAssignments)
}
## =============================================================================

update_cluster_labels <- function(oldSeqsClustLabels, collatedClustAssignments,
                                  flags) {
    # print("=== updating clust labels ===")
    # print(collatedClustAssignments)
    # print("=== oldSeqLustLabels ===")
    # print(oldSeqsClustLabels)
    nClusters <- length(collatedClustAssignments)
    # print(paste0("nClusters: ", nClusters))
    if (flags$verboseFlag) {
        cat("Updating sequence cluster labels\n")
        cat(paste0("#Clusters: ", length(collatedClustAssignments), "\n"))
    }
    newSeqsClustLabels <- oldSeqsClustLabels
    if (nClusters > 9) {
        message("More than 9 clusters")
    }
    for (i in 1:length(collatedClustAssignments)) {
        needUpdateIdx <- collatedClustAssignments[[i]]
        newSeqsClustLabels[needUpdateIdx] <- sapply(
            newSeqsClustLabels[needUpdateIdx],
            function(x) {
                paste0(c(x, toString(i)), collapse = "-")
            })
    }
    # print("=== newSeqClustLabels ===")
    # print(newSeqsClustLabels)
    return(newSeqsClustLabels)
}
## =============================================================================

prepare_chunks <- function(total_avail, reqdChunkSize, checkLength,
                           flags = list(debugFlag = TRUE,
    verboseFlag = TRUE, plotVerboseFlag = TRUE, timeFlag = TRUE)) {
    ### total_avail is the set of seq_ids to be chunked (not array indices)
    if (length(total_avail) > reqdChunkSize) {
        chunkStarts <- seq(1, length(total_avail), by = reqdChunkSize)
        chunkEnds <- seq(reqdChunkSize, length(total_avail), by = reqdChunkSize)
        if (length(chunkStarts) > length(chunkEnds)) {
            if (flags$debugFlag) {
                print("Chunk starts/ends altered")
            }
            chunkStarts <- chunkStarts[-c(length(chunkStarts))]
            chunkEnds[length(chunkEnds)] <- checkLength
        }
        if (flags$debugFlag) {
            print(chunkStarts)
            print(chunkEnds)
        }
        preparedChunks <- vector("list", length(chunkStarts))
        for (i in 1:length(chunkStarts)) {
            preparedChunks[[i]] <- total_avail[chunkStarts[i]:chunkEnds[i]]
        }
    } else {
        preparedChunks <- vector("list", 1)
        preparedChunks[[1]] <- total_avail
    }
    return(preparedChunks)

}
## =============================================================================

handle_clustering_of_factors <- function(globFactorsMat,
                                         distMethod = "cosangle",
                                         flags) {
    ### Currently relying on HOPACH algorithm - Compute cosine similarities
    ### (values between 0-1) - Using HOPACH algorithm for clustering with chosen
    ###  distance measure (we currently use 'cosangle' distance measure)
    ###  globFactorsDistMat <- hopach::distancematrix(t(globFactorsMat),
    ###  d = distMethod)
    globFactorsDistMat <- compute_factor_distances(globFactorsMat,
                                                   distMethod = distMethod)

    doHopach <- hopach::msscheck(globFactorsDistMat, within = "mean",
                                  between = "mean")

    if (doHopach[1] == 1 && ncol(globFactorsMat) > 2) {
        message("Error expected to occur now")
    }
    #
    globFactorsHopach <- hopach::hopach(data = t(globFactorsMat),
                                        dmat = globFactorsDistMat,
                                        d = "cosangle", newmed = "uwnn",
                                        clusters = "best", coll = "all",
                                        initord = "clust", mss = "mean",
                                        verbose = flags$verboseFlag)
    if (flags$debugFlag || flags$verboseFlag) {
        cat(paste0("HOPACH identified #clusters:",
                   globFactorsHopach$clustering$k,
            "\n"))
    }
    if (flags$plotVerboseFlag) {
        ### Order the medians, accordingly change order of the collated cluster
        ### assignments
        medoidsIdx <- get_hopach_cluster_medoidsIdx(globFactorsHopach)
        print(medoidsIdx)
    }
    return(globFactorsHopach)
}
## =============================================================================

# final_adjustment_of_ordering <- function(archRresult) {
#     # Finally, arrange clusters from processed outer chunks.
#
#
#
# }
## =============================================================================
