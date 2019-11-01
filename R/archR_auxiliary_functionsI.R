#' @title Get cluster labels for chosen level/iteration
#'
#' @description Given a seqsClustLabels, collect cluster labels for sequences at
#'  the chosen iteration/level.
#'
#' @param given_seqsClustLabels from archR result object
#' @param chooseLevel choose a level/iteration. This value is number of
#' iterations + 1
#'
#' @return A numeric vector of the same size as seqsClustLabels, with labels
#' only up to the chosen iteration
#' @export
#'
collect_cluster_labels <- function(given_seqsClustLabels, chooseLevel = 1) {
    ## Check if all_ok, all elements should have same length
    .assert_archR_seqsClustLabels_at_end(given_seqsClustLabels)
    splitChar <- "-"
    elements_length <- unique(unlist(lapply(strsplit(given_seqsClustLabels,
                                                    split = splitChar),
                                                    length)))
    if (chooseLevel > elements_length) {
        stop("Choose_levels (", chooseLevel,
        ") greater than levels present in cluster_labels(",
        elements_length, ").")
    } else {
        selectedLabels <- unlist(lapply(strsplit(given_seqsClustLabels,
                                                split = splitChar),
                                        function(x) {
                                            paste0(x[seq_len(chooseLevel)],
                                                collapse = splitChar)
                                        }))
    }
    .assert_archR_seqsClustLabels_at_end(selectedLabels)
    return(selectedLabels)
}
## =============================================================================

# @title Assign samples to clusters
#
# @param samplesMatrix
#
# @return Similar to .map_clusters_to_factors, A list that can be assigned as
# an element in globClustAssignments
#
.assign_samples_to_clusters <- function(samplesMatrix, iChunksColl,
                                        iChunkIdx) {
    .assert_archR_samplesMatrix(samplesMatrix)
    nClusters = nrow(samplesMatrix)
    returnClusterAsList <- vector("list", nClusters)
    ## Fetch the cluster memberships for each sample (along the columns)
    clustMemberships <- apply(samplesMatrix, 2, which.max)
    if (length(unique(clustMemberships)) != nClusters) {
        print(unique(clustMemberships))
        print(length(unique(clustMemberships)))
        print(nClusters)
        stop("Basis vector that got no sequences assigned: ",
             setdiff( unique(clustMemberships), seq_len(nClusters)))
    } else {
        for (i in seq_along(returnClusterAsList)) {
            returnClusterAsList[[i]] <-
                iChunksColl[[iChunkIdx]][clustMemberships == i]
        }
        return(returnClusterAsList)
    }
}


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
            warning(c("Factor(s) got no sequences assigned: ",
                    paste0(thisGotLeftOut, sep="-")), immediate. = TRUE)
        }
        return(rightClusterOrders)
    }
}
## =============================================================================


## @title Collate clusters
##
## @param hopachObj A hopach object which is a special list with elements
## having fixed names
## @param globClustAssignments A list
##
## @return collatedClustAssignments A list
.collate_clusters <- function(hopachObj, globClustAssignments) {
    ## Based on which factors are being combined (hopachObj), look at the
    ## globClustAssignments variable to combine the respective sequences
    ## together.
    .assert_archR_globClustAssignments(globClustAssignments)
    if (is.null(hopachObj)) {
        ## When HOPACH clustering was not performed, the globClustAssingments
        ## variable is directly assigned
        collatedClustAssignments <- globClustAssignments
    } else {
        ## When HOPACH clustering was performed, liase with the hopachObj
        .assert_archR_hopachObj(hopachObj, test_null = FALSE)
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
                            function(x) {globClustAssignments[[x]]}
                            )
                    )
            }
        }  ## collation for loop ends
    }
    return(collatedClustAssignments)
}
## =============================================================================

## @title Update labels of sequences in a cluster
## @param oldSeqsClustLabels
##
## @param collatedClustAssignments
## @param flags List. Flags variable from archR config
##
## @return newSeqsClustLabels Vector os clustLabels
.update_cluster_labels <- function(oldSeqsClustLabels, collatedClustAssignments,
                                    flags = list(debugFlag = FALSE,
                                                verboseFlag = TRUE,
                                                plotVerboseFlag = FALSE,
                                                timeFlag = FALSE)) {
    .assert_archR_seqsClustLabels(oldSeqsClustLabels)
    .assert_archR_globClustAssignments(collatedClustAssignments)
    .assert_archR_flags(flags)
    # .assert_archR_consistent_nSeqs_w_clusters(oldSeqsClustLabels,
    #                                     collatedClustAssignments)
    nClusters <- length(collatedClustAssignments)
    ## When #clusters > 9, we use LETTERS as cluster labels instead of numbers
    ## 10, 11, 12 and so. -- This does not work when nClusters > 50 or 100 or
    ## so.
    ## Instead, use numerics w/ as.character and pre-sort to have them in the
    ## order that will be returned by levels (in get_seqs_clusters_in_list fn)
    candidateClustLabels <- sort(as.character(seq_len(nClusters)))
    if (flags$verboseFlag) {
        message("Updating sequence cluster labels")
        message("#Clusters: ", nClusters)
    }
    ## This was needed when cluster labels had no separator. May be removed now.
    ## if (nClusters > 9) {
    ##     warning("More than 9 clusters")
    ## }
    newSeqsClustLabels <- oldSeqsClustLabels
    for (i in seq_len(nClusters)) {
        needUpdateIdx <- collatedClustAssignments[[i]]
        newSeqsClustLabels[needUpdateIdx] <-
            vapply(newSeqsClustLabels[needUpdateIdx],
                    function(x) {
                        paste0(c(x, candidateClustLabels[i]), collapse = "-")
                    }, character(1)
                )
    }
    .assert_archR_seqsClustLabels(newSeqsClustLabels)
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
                                                    verboseFlag = FALSE,
                                                    plotVerboseFlag = FALSE,
                                                    timeFlag = FALSE)
                                        ) {
    .assert_archR_featuresMatrix(globFactorsMat)
    .assert_archR_flags(flags)
    ##
    globFactorsDistMat <- .compute_factor_distances(globFactorsMat,
                                                distMethod = distMethod)
    ## These lines were put for additional caution/check. May not be needed now.
    ## doHopach <- hopach::msscheck(globFactorsDistMat, within = "mean",
    ##                             between = "mean")
    ## if (doHopach[1] == 1 && ncol(globFactorsMat) > 2) {
    ##     message("Error expected to occur now")
    ## }
    ##
    globFactorsHopach <- hopach::hopach(data = t(globFactorsMat),
                                        dmat = globFactorsDistMat,
                                        d = "cosangle", newmed = "uwnn",
                                        clusters = "best", coll = "all",
                                        initord = "clust", mss = "mean",
                                        verbose = flags$verboseFlag)
    if (flags$debugFlag ||
        flags$verboseFlag) {
        message("Identified #clusters: ", globFactorsHopach$clustering$k)
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

#' Keep this internal or external function?
#' @description We use hclust/hierarchical clustering for reordering archR
#' clusters
#'
#' @title Reorder archR clustering at the current last level
#'
#' @param archRresult The archRresult object
#'
#' @return Ordering returned from hclust
#'
#' @importFrom stats hclust dist
#' @export
reorder_archRresult <- function(archRresult) {
    # Depends on archRresult object having a fixed set of names.
    # We need to .assert them
    # Finally, arrange clusters from processed outer chunks.
    #  using hclust
    lastLevel <- length(archRresult$clustBasisVectors)
    temp_hclust <-
        stats::hclust(stats::dist(
            t(archRresult$clustBasisVectors[[lastLevel]]$basisVectors)),
                             method = "ave")
    new_order <- temp_hclust$order
    # Fetch original seqClustLabels as a list
    origSeqsClustersAsList <- get_seqs_clusters_in_a_list(
                        archRresult$seqsClustLabels,
                        chooseLevel = lastLevel + 1)
    if (length(new_order) != length(origSeqsClustersAsList)) {
        print("===== SAMARTH SAMARTH =====")
        print(archRresult$seqsClustLabels)
        print(origSeqsClustersAsList)
        print(length(origSeqsClustersAsList))
        print(new_order)
        print(length(new_order))
        stop("SAMARTH: Error")
    }
    ## arrange by the new ordering
    newSeqsClusters <- lapply(new_order, function(x){
                                            origSeqsClustersAsList[[x]]
                                            }
                            )
    ## the labels should be set as.character
    newSeqsClustLabels <- unlist(lapply(
                                    seq_along(newSeqsClusters),
                                    function(x){
                                        as.character(
                                        rep(x, length(newSeqsClusters[[x]]))
                                        )
                                    }
                                )
                            )
    ##
    newClustBasisVectors <-
        archRresult$clustBasisVectors[[lastLevel]]$basisVectors[, new_order]
    ##
    new_field <- list(seqsClusters = newSeqsClusters,
                        seqsClustLabels = newSeqsClustLabels,
                        clustBasisVectors = newClustBasisVectors)
    ##
    archRresult$final <- new_field
    ##
    return(archRresult)
}
## =============================================================================

#' @title Retrieve sequence clusters as a list
#'
#' @param seqsClustLabels Sequences with cluster labels as in the archR result
#' object
#' @param chooseLevel Specify the level (archR iteration) at which sequence
#' clusters are to be reported. Default is 1.
#'
#' @return A list holding sequence clusters
#' @export
get_seqs_clusters_in_a_list <- function(seqsClustLabels, chooseLevel = 1){

    chosenLevelLabels <- collect_cluster_labels(seqsClustLabels,
                                                chooseLevel = chooseLevel)
    clusterLevels <- levels(as.factor(chosenLevelLabels))
    # message(length(clusterLevels),
    #         " clusters identified by labels: ",
    #         paste0(clusterLevels, sep = "**"))
    seqs_clusters_as_a_list <- sapply(seq_along(clusterLevels),
                                        function(x){
                                            which(chosenLevelLabels ==
                                                    clusterLevels[x])
                                    })

    return(seqs_clusters_as_a_list)
}
## =============================================================================
