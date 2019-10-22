#' @title Get cluster labels for chosen level/iteration
#'
#' @description Given a seqsClustLabels, collect cluster labels for sequences at
#'  the chosen iteration/level.
#'
#' @param given_seqsClustLabels from archR result object
#' @param choose_levels choose a level/iteration
#'
#' @return A numeric vector of the same size as seqsClustLabels, with labels
#' only up to the chosen iteration
#' @export
#'
collect_cluster_labels <- function(given_seqsClustLabels, choose_levels = 1) {
    ## Check if all_ok, all elements should have same length
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
                paste0(x[seq_len(choose_levels)], collapse = splitChar)
            }))
    }
    return(selectedLabels)
}
## =============================================================================

# @title Map clusters to factors
#
# @description This function maps the clusters with the corresponding factors
#  factors <--> clusters; meaning if there are 5 clusters and 5 factors, which
#  cluster corresponds to which factor is assigned by this function.
#  It takes as input the clusters, and returns the correct order of the
#  clusters.
#
# @param samplesMatrix A matrix. Samples matrix from solving NMF
# @param clustOrderIdx A nested list similar to reordering_idx from kmeans
# @param iChunksColl A list
# @param iChunkIdx A numeric
# @param flags List with four fixed elements. flags param in archR config
#
# @return A list that can be assigned as an element in globClustAssignments
.map_clusters_to_factors <- function(samplesMatrix, clustOrderIdx, iChunksColl,
                                    iChunkIdx, flags) {
    ## clustOrderIdx IS clustering_sol_kmeans$reordering_idx, the reorderingIdx
    ## is a nested list, i.e. a list holding list of seqs_ids belonging to one
    ## cluster together. Hence, length(this variable) gives #clusters.
    ## samplesMatrix IS samplesMatrix.
    ##
    ## Make checks
    .assert_archR_samplesMatrix(samplesMatrix)
    ##
    if (!is.list(clustOrderIdx)) {
        stop("Mapping clusters to factors, cluster orders is not a list")
    } else if (!is.vector(clustOrderIdx[[1]])) {
            stop("Mapping clusters to factors, cluster orders is not a vector")
    }
    ## check iChunksColl and iChunkIdx
    if (!is.null(iChunksColl) && length(iChunksColl) < 1) {
        stop("Mapping clusters to factors, expecting inner chunks as a
                list of length > 0")
    }
    if (!is.null(iChunkIdx) && !is.numeric(iChunkIdx) && iChunkIdx > 0) {
        stop("Mapping clusters to factors, expecting inner chunk index as a
             numeric > 0")
    }
    ##
    .assert_archR_flags(flags)
    ##
    rightClusterOrders <- vector("list", length(clustOrderIdx))
    if (length(clustOrderIdx) == 1) {
        ## Special case
        rightClusterOrders[[1]] <- iChunksColl[[iChunkIdx]][clustOrderIdx[[1]]]
        return(rightClusterOrders)
        ##
    } else if (length(clustOrderIdx) > 1) {
        ##
        for (cluster_idx in seq_along(clustOrderIdx)) {
            relevant_factor <-
                which.max(rowMeans(
                as.matrix(samplesMatrix[, clustOrderIdx[[cluster_idx]]])
                ))
            if (length(relevant_factor) > 1) {
                if (flags$debugFlag) {
                    message("Multiple basis vectors attained max.")
                    message(relevant_factor)
                }
            } else {
                if (flags$debugFlag) {
                    message("Factor-mapping OK")
                    message(relevant_factor)
                }
                rightClusterOrders[[relevant_factor]] <-
                    iChunksColl[[iChunkIdx]][clustOrderIdx[[cluster_idx]]]
            }
        }
        ## Check if any factor got left out?
        ## i.e. no seqs assigned to its cluster
        if (any(lapply(rightClusterOrders, length) == 0)) {
            thisGotLeftOut <- which(lapply(rightClusterOrders, length) == 0)
            message("WARNING: Factor(s) got no sequences assigned: ",
                    thisGotLeftOut)
            warning("WARNING: Factor(s) got no sequences assigned: ",
                    thisGotLeftOut)
        }
        return(rightClusterOrders)
    }
}
## =============================================================================

.collate_clusters <- function(hopachObj, globClustAssignments) {
    ## Based on which factors are being combined (hopachObj), look at the
    ## globClustAssignments variable to combine the respective sequences
    ## together.
    if (is.null(hopachObj)) {
        # When HOPACH not performed, the global_cluster_assingments variable
        # holds a list of lists, so we should flatten it with unlist()
        collatedClustAssignments <- globClustAssignments
        ##
    } else {
        ##
        nClusters <- hopachObj$clustering$k
        clustSizes <- hopachObj$clustering$sizes
        elementsOrder <- hopachObj$clustering$order
        ##
        collatedClustAssignments <- vector("list", nClusters)
        ##
        for (coll_cluster_idx in seq_along(collatedClustAssignments)) {
            to_combine_this_itr <- clustSizes[coll_cluster_idx]
            leave_first <- 0
            if (coll_cluster_idx > 1) {
                leave_first <- sum(clustSizes[seq_len(coll_cluster_idx - 1)])
            }
            pick_this_itr <- elementsOrder[
                (leave_first + 1):(leave_first + to_combine_this_itr)]
            ##
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
        }  ## collation for loop ends
    }
    return(collatedClustAssignments)
}
## =============================================================================

.update_cluster_labels <- function(oldSeqsClustLabels, collatedClustAssignments,
                                    flags) {
    nClusters <- length(collatedClustAssignments)
    if (flags$verboseFlag) {
        cat("Updating sequence cluster labels\n")
        cat(paste0("#Clusters: ", length(collatedClustAssignments), "\n"))
    }
    newSeqsClustLabels <- oldSeqsClustLabels
    if (nClusters > 9) {
        message("More than 9 clusters")
    }
    for (i in seq_along(collatedClustAssignments)) {
        needUpdateIdx <- collatedClustAssignments[[i]]
        newSeqsClustLabels[needUpdateIdx] <- vapply(
            newSeqsClustLabels[needUpdateIdx],
            function(x) {
                paste0(c(x, toString(i)), collapse = "-")
            }, character(1))
    }
    return(newSeqsClustLabels)
}
## =============================================================================

# @title Prepare chunks out of given set of sequences (sequence IDs)
# @param total_set set of sequences to be chunked
# @param reqdChunkSize the given chunk size
# @param flags Specify the flags from the config param set for archR
#
# @return A list with chunks of sequences (sequence IDs) as its elements
.prepare_chunks <- function(total_set, reqdChunkSize,
                            flags = list(debugFlag = FALSE,
                                        verboseFlag = TRUE,
                                        plotVerboseFlag = FALSE,
                                        timeFlag = FALSE)
                            ) {
    ## total_set is the set of seq_ids to be chunked (not array indices)
    if (is.null(total_set)) {
        stop("Preparing chunks, 'total_set' is NULL")
    }
    chunkLength <- length(total_set)
    if (chunkLength == 0) {
        stop("Preparing chunks, length of 'total_set' is 0")
    }
    .assert_archR_innerChunkSize_independent(reqdChunkSize)
    .assert_archR_flags(flags)
    ##
    if (chunkLength > reqdChunkSize) {
        chunkStarts <- seq(1, chunkLength, by = reqdChunkSize)
        chunkEnds <- seq(reqdChunkSize, chunkLength, by = reqdChunkSize)
        if (length(chunkStarts) > length(chunkEnds)) {
            if ((chunkLength - chunkStarts[length(chunkStarts)]) >
                round(0.5*reqdChunkSize)) {
                chunkEnds <- append(chunkEnds, chunkLength)
            } else {
                chunkStarts <- chunkStarts[-length(chunkStarts)]
                chunkEnds[length(chunkEnds)] <- chunkLength
            }
        }
        if (flags$debugFlag) {
            print(chunkStarts)
            print(chunkEnds)
        }
        preparedChunks <- vector("list", length(chunkStarts))
        for (i in seq_along(chunkStarts)) {
            preparedChunks[[i]] <- total_set[chunkStarts[i]:chunkEnds[i]]
        }
    } else {
        preparedChunks <- vector("list", 1)
        preparedChunks[[1]] <- total_set
    }
    if (length(preparedChunks) < 1) {
        stop("Preparing chunks, length was 0")
    }
    return(preparedChunks)

}
## =============================================================================

# @title Handle clustering of NMF factors
#
# @param globFactorsMat Specify the NMF factors as a matrix with the individual
# factors along the columns.
# @param distMethod Specify the distance measure to be used. Default is cosine
# measure.
# @param flags Specify the flags from the config param set for archR
#
# @return hopach object
.handle_clustering_of_factors <- function(globFactorsMat,
                                        distMethod = "cosangle",
                                        flags = list(debugFlag = FALSE,
                                                    verboseFlag = TRUE,
                                                    plotVerboseFlag = FALSE,
                                                    timeFlag = FALSE)
                                        ) {
    .assert_archR_flags(flags)
    ## Currently relying on HOPACH algorithm - Compute cosine similarities
    ## (values between 0-1) - Using HOPACH algorithm for clustering with chosen
    ##  distance measure (we currently use 'cosangle' distance measure)
    globFactorsDistMat <- .compute_factor_distances(globFactorsMat,
                                                distMethod = distMethod)
    ##
    doHopach <- hopach::msscheck(globFactorsDistMat, within = "mean",
                                between = "mean")

    if (doHopach[1] == 1 && ncol(globFactorsMat) > 2) {
        message("Error expected to occur now")
    }
    ##
    globFactorsHopach <- hopach::hopach(data = t(globFactorsMat),
                                        dmat = globFactorsDistMat,
                                        d = "cosangle", newmed = "uwnn",
                                        clusters = "best", coll = "all",
                                        initord = "clust", mss = "mean",
                                        verbose = flags$verboseFlag)
    if (flags$debugFlag ||
        flags$verboseFlag) {
        cat(paste0(
            "Identified #clusters:",
            globFactorsHopach$clustering$k,
            "\n"
        ))
    }
    if (flags$plotVerboseFlag) {
        ## Order the medians, accordingly change order of the collated cluster
        ## assignments
        medoidsIdx <-
            .get_hopach_cluster_medoidsIdx(globFactorsHopach)
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
