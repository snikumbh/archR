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
.get_cluster_memberships_per_run <- function(samplesMatrix, iChunksColl,
                                        iChunkIdx) {
    .assert_archR_samplesMatrix(samplesMatrix)
    nClusters = nrow(samplesMatrix)
    returnClusterAsList <- vector("list", nClusters)
    ## Fetch the cluster memberships for each sample (along the columns)
    clustMemberships <- apply(samplesMatrix, 2, which.max)
    ##
    if (length(unique(clustMemberships)) != nClusters) {
        # print(unique(clustMemberships))
        # print(length(unique(clustMemberships)))
        # print(nClusters)
        message("WARNING: Basis vector that got no sequences assigned: ",
             setdiff( unique(clustMemberships), seq_len(nClusters)))
    } else {
        ##
    }
    return(clustMemberships)
}


.assign_samples_to_clusters <- function(clusterMembershipsVec, nClusters,
                                        iChunksColl, iChunkIdx) {
    returnClusterAsList <- vector("list", nClusters)
    for (i in seq_along(returnClusterAsList)) {
        returnClusterAsList[[i]] <-
            iChunksColl[[iChunkIdx]][clusterMembershipsVec == i]
    }
    return(returnClusterAsList)
}



.get_consensusClustMemberships <- function(allRunMemberships, nClusters) {

    consensusClustMemberships <-
        apply(allRunMemberships, 2, function(yVec){
                        which.max(
                            vapply(seq_len(nClusters),
                                   function(x){sum(yVec == x)},
                            numeric(1))
                        )}
        )

    return(consensusClustMemberships)
}



.reorderFeatureMatsInList <- function(featureMatsList) {
    ## Look at the first matrix and order other according to it
    forOrderCorrespondance <- as.matrix(featureMatsList[[1]])
    reorderedFeatureMatsList <- vector("list", length(featureMatsList))
    reorderingIdx <- vector("list", length(featureMatsList))
    ## use cosine distance to map them
    for (ll in seq_along(featureMatsList)) {
        # for (i in 1:ncol(forOrderCorrespondance)) {
        #     temp <- hopach::distancevector(t(featureMatsList[[ll]]),
        #                     as.vector(forOrderCorrespondance[,i]),
        #                     d = "cosangle")
        #     print(temp)
        # }
        # print("NEXT RUN")
        factorOrder <- apply(forOrderCorrespondance, 2, function(x){
                                which.min(
                                    hopach::distancevector(
                                    t(featureMatsList[[ll]]),
                                    as.vector(x), d = "cosangle")
                                )
                                }
                            )
        # print(factorOrder)
        # print(ncol(forOrderCorrespondance))
        if (length(unique(factorOrder)) != ncol(forOrderCorrespondance)) {
            print(factorOrder)
            print(ncol(forOrderCorrespondance))
        }
        #stopifnot(length(unique(factorOrder)) == ncol(forOrderCorrespondance))
        reorderedFeatureMatsList[[ll]] <- featureMatsList[[ll]][ , factorOrder]
        reorderingIdx[[ll]] <- factorOrder
    }
    #stop("SAMARTH")
    return(list(reorderedFeatureMatsList = reorderedFeatureMatsList,
                reorderingIdx = reorderingIdx))
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
        print(clustOrderIdx)
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
        message("Cluster to Factor:")
        for (cluster_idx in seq_along(clustOrderIdx)) {
            ##
            # print(apply(as.matrix(samplesMatrix[, clustOrderIdx[[cluster_idx]]]),
            #             1, quantile))
            # print(as.matrix(samplesMatrix[, clustOrderIdx[[cluster_idx]]]))
            rowwiseMeans <- rowMeans(
                as.matrix(samplesMatrix[, clustOrderIdx[[cluster_idx]]])
            )
            print(paste0("Cluster ID:", cluster_idx))
            print(paste0("Rowwise means", rowwiseMeans))
            relevant_factor <- which.max(rowwiseMeans)
            relevant_descending_order <-
                sort(samplesMatrix[relevant_factor,clustOrderIdx[[cluster_idx]]],
                     index.return = TRUE,
                     decreasing = TRUE)$ix
            # print(relevant_factor)
            # #stop("SAMARTH")
            ##
            # relevant_factor <-
            #     which.max(rowMeans(
            #     as.matrix(samplesMatrix[, clustOrderIdx[[cluster_idx]]])
            #     ))
            if (length(relevant_factor) > 1) {
                if (flags$debugFlag) {
                    message("Multiple basis vectors attained max.")
                    message(relevant_factor)
                }
            } else {
                message(relevant_factor)
                if (flags$debugFlag) {
                    message("Factor-mapping OK")
                    message(relevant_factor)
                }

                seqs_vector_to_assign <-
                    iChunksColl[[iChunkIdx]][clustOrderIdx[[cluster_idx]]]
                ## order the sequences by the scores (highest first, weakest last)
                rightClusterOrders[[relevant_factor]] <-
                    seqs_vector_to_assign[relevant_descending_order]
            }
        }
        ## Check if any factor got left out?
        ## i.e. no seqs assigned to its cluster
        if (any(lapply(rightClusterOrders, length) == 0)) {
            print(clustOrderIdx)
            thisGotLeftOut <- which(lapply(rightClusterOrders, length) == 0)
            warning(c("Factor(s) got no sequences assigned: ",
                    paste0(thisGotLeftOut, sep="-")), immediate. = TRUE)
            for (id in thisGotLeftOut) {
                rightClusterOrders[[id]] <- 0
            }
        }
        return(rightClusterOrders)
    }
}
## =============================================================================

## @title Collate clusters from hierarchical clustering
##
## @param clusterings A list object giving indices of factors in a clusters
## @param globClustAssignments A list
##
## @return collatedClustAssignments A list
.collate_clusters2 <- function(clustList, globClustAssignments,
                               flags = list(debugFlag = FALSE,
                                            verboseFlag = FALSE,
                                            plotVerboseFlag = FALSE,
                                            timeFlag = FALSE)) {
    # print(globClustAssignments)
    if(flags$debugFlag){
        message("=== ^ That is what I got ===")
    }
    .assert_archR_globClustAssignments(globClustAssignments)
    nClusters <- length(clustList)
    clustSizes <- unlist(lapply(clustList, length))
    ##
    collatedClustAssignments <- vector("list", nClusters)
    for(clustIdx in seq_along(collatedClustAssignments)){
        temp <- unlist(lapply(clustList[[clustIdx]], function(x){
            globClustAssignments[[x]]
        }))
        if(flags$debugFlag){
            message("=== INFO ===")
            print(clustList[[clustIdx]])
            message("Size:", length(temp))
        }
        collatedClustAssignments[[clustIdx]] <- temp
    }
    if(flags$debugFlag){
        message("=== I am returning ===")
        print(collatedClustAssignments)
    }

    return(collatedClustAssignments)
}



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
    # print("Candidate clust labels")
    # print(candidateClustLabels)
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
        # print("Need update idx")
        # print(head(needUpdateIdx))
        # print(length(needUpdateIdx))
        newSeqsClustLabels[needUpdateIdx] <-
            vapply(newSeqsClustLabels[needUpdateIdx],
                    function(x) {
                        #paste0(c(x, candidateClustLabels[i]), collapse = "-")
                        paste0(candidateClustLabels[i])
                    }, character(1)
                )
    }
    ##
    .assert_archR_seqsClustLabels(newSeqsClustLabels)
    # print("Inside update cluster labels")
    # print(oldSeqsClustLabels)
    # print(newSeqsClustLabels)
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
                                        clustMethod = "hc",
                                        linkage = "average",
                                        distMethod = "cosangle",
                                        flags = list(debugFlag = FALSE,
                                                    verboseFlag = FALSE,
                                                    plotVerboseFlag = FALSE,
                                                    timeFlag = FALSE),
                                        returnOrder = TRUE) {
    .assert_archR_featuresMatrix(globFactorsMat)
    .assert_archR_flags(flags)
    ##
    globFactorsDistMat <- .compute_factor_distances(globFactorsMat,
                                                distMethod = distMethod)

    if(clustMethod == "hopach"){
        ## These lines were put for additional caution/check.
        ## May not be needed now.
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
    } else if(clustMethod == "hc"){
        if(flags$debugFlag) message("Hierarchical clustering instead of HOPACH")
        ############### REORDER CLUSTERS FROM OUTER CHUNK ######################
        if(flags$debugFlag) message("SAMARTH: RE-ORDERING CLUSTERS")

        if(flags$verboseFlag) {
            message("Cluster order by ", linkage," linkage w/ ",
                    distMethod, " distance:")
        }
        ##
        as_dist_mat <- stats::as.dist(hopach::as.matrix(globFactorsDistMat))
        temp_hclust <- stats::hclust(as_dist_mat, method = linkage)
        new_order <- temp_hclust$order
        ##
        if(flags$debugFlag){
            print(new_order)
            for(k in 1:(length(new_order)-1) ){
                message(new_order[k], ", ", new_order[k+1],", ",
                        globFactorsDistMat[new_order[k], new_order[k+1]])
            }
            print(new_order)
            ##
            message("=== CALLING FUNCTION === ")
        }
        ##

        if(returnOrder){
            ## return just the new ordering
            return(new_order)
        }else{
            clustList <- .get_clusters_from_hc(hcObj = temp_hclust,
                                               distMat = globFactorsDistMat,
                                               verbose = flags$debugFlag)
            if(flags$debugFlag) {
                message("=== DONE ===")
                print(clustList)
            }
            return(clustList)
        }
    }
}
## =============================================================================

# .get_clusters_from_hc2 <- function(hcObj, distMat, verbose = FALSE){
#
#
# }



.get_clusters_from_hc <- function(hcObj, distMat, verbose = FALSE){

    #####
    ## Algo to cut dendrogram and form clusters
    ## corner case: 1. what if start nodes are not available?
    if(verbose){
        message("=== *** Merge ===")
        print(hcObj$merge)
        message("=== *** height ===")
        print(hcObj$height)
        message("=== *** dendro ===")
        print(utils::str(stats::as.dendrogram(hcObj)))
        message("=== *** ===")
    }
    newOrder <- hcObj$order

    leftN <- rep(NA, length(newOrder))
    rightN <- rep(NA, length(newOrder))
    chooseVal <- rep(NA, length(newOrder))
    if(verbose) {
        message("Length: ", length(newOrder))
    }
    for(k in seq_along(newOrder)){
        idxInOrder <- which(newOrder == k)

        if(idxInOrder != 1) leftN[k] <- newOrder[idxInOrder - 1]
        if(idxInOrder != length(newOrder)) rightN[k] <- newOrder[idxInOrder + 1]
        if(!is.na(leftN[k]) && !is.na(rightN[k]) ){
            if(distMat[leftN[k],k] > distMat[k,rightN[k]] ){
                chooseVal[k] <- rightN[k]}
            else {
                chooseVal[k] <- leftN[k]}
        }
        #
        specialIdx <- which(is.na(leftN))
        chooseVal[specialIdx] <- rightN[specialIdx]
        specialIdx <- which(is.na(rightN))
        chooseVal[specialIdx] <- leftN[specialIdx]
        #
        if(verbose){
        message("[", leftN[k], ",", rightN[k], "]",
                " for ", k, ", chose ", chooseVal[k])
        }
    }

    # find start nodes
    startNodeList <- list()
    chosenIdx <- vector() #stores chooseVal indices to ignore later
    for(k in 1:length(chooseVal)){
        if(chooseVal[chooseVal[k]] == k){
            chosenIdx <- append(chosenIdx, k)
            if(k < chooseVal[k]){
                a <- k
                b <- chooseVal[k]
            } else{
                a <- chooseVal[k]
                b <- k
            }
            nodeEntry <- c(a,b)
            matches <- unlist(lapply(startNodeList, function(l) {nodeEntry == l}))
            if(!(any(matches))){
                startNodeList <- append(startNodeList, NA)
                startNodeList[[length(startNodeList)]] <- nodeEntry
                if(verbose){
                    message(a,"-",b)
                }
            }
        }
    }
    relIdx <- setdiff(1:length(chooseVal), chosenIdx)

    # case when all nodes end up as part of startNodes,
    # relidx is numeric(0), length 0
    nodeList <- startNodeList
    if(length(relIdx) != 0 && length(relIdx) < length(chooseVal)){
        ##
        reChooseVal <- chooseVal
        reChooseVal[-relIdx] <- NA
        if(verbose){
            message("Stringing nodes on the left")
            print(reChooseVal)
            message("=== Start nodes: ===")
            print(startNodeList)
            message("=== Done ===")
        }
        for(l in 1:length(startNodeList)){
            # if(verbose) message(nodeList[[l]])
            a <- startNodeList[[l]][1]
            # if(verbose) message(a, ", ", b)
            while(any(which(reChooseVal == a)) ){
                # if (verbose) print(reChooseVal)
                matchIdx <- which(reChooseVal == a)
                nodeList[[l]] <- append(nodeList[[l]], matchIdx, after = 0)
                a <- matchIdx
                # if (verbose) print(a)
                reChooseVal[matchIdx] <- NA
                # if (verbose) print(reChooseVal)
            }
        }
        if(verbose){
            # print(reChooseVal)
            # message("=== Nodes list: ===")
            # print(nodeList)
            # message("=== Done ===")
        }
        if(!all(is.na(reChooseVal))){
            if(verbose){
                message("Stringing nodes on the right")
            }
            for(l in 1:length(nodeList)){
                if(verbose) message(paste0(nodeList[[l]], sep=","))
                b <- startNodeList[[l]][2]
                # if(verbose) message(a, ", ", b)
                while(any(which(reChooseVal == b)) ){
                    if (verbose) print(reChooseVal)
                    matchIdx <- which(reChooseVal == b)
                    nodeList[[l]] <- append(nodeList[[l]], matchIdx,
                                            after = length(nodeList[[l]]))
                    b <- matchIdx
                    if (verbose) print(b)
                    reChooseVal[matchIdx] <- NA
                    if (verbose) print(reChooseVal)
                }
            }
        }
    }
    if (verbose){
        print(nodeList)
        message("Re-adjust order")
    }
    starts <- lapply(nodeList, function(x){utils::head(x, 1)})
    matchOrder <- vapply(starts, function(x){which(x == newOrder)}, numeric(1))
    newIdx <- sort(matchOrder, index.return = TRUE)$ix
    nodeListUpd <- lapply(newIdx, function(x){nodeList[[x]]})
    nodeList <- nodeListUpd
    if (verbose) print(nodeList)

    if(verbose){
        message("Final Length: ", sum(unlist(lapply(nodeList, length))))
    }
    stopifnot(sum(unlist(lapply(nodeList, length))) == length(newOrder))
    #####
    return(nodeList)
}


# compare_NMF_factors <- function(givenFactors){
#
#     archRresult_samarth <- readRDS(file.path("/mnt/biggley/home/sarvesh/projects/promoter-architectures-all/clustering-promoter-architectures/promArch/archR/vignettes/archR_on_dinucleotide_profiles",
#                                      "experiments-hierarchical-for-hopach/run-one_43/archRresult.rds"))
#
#
#     sam_mat <- archRresult_samarth$clustBasisVectors[[3]]$basisVectors
#
#     dim_names <- archR::get_dimers_from_alphabet(c("A", "C", "G", "T"))
#     sam_mat2 <- lapply(1:ncol(sam_mat),
#                        function(x){matrix(sam_mat[,x], nrow = 1456/91, byrow = TRUE,
#                                           dimnames = list(dim_names))
#                            })
#     sam_mat2 <- lapply(sam_mat2, function(x){
#         colnames(x) <- positions
#         x})
#     ## Go over windows of a small size (say, 5 bp)
#     ## What are the top-3 dimers in this window?
#     window_size <- 5
#     windowStarts <- seq(1, 91-window_size+1, by = 1)
#     windowEnds <- seq(window_size, 91, by = 1)
#
#     # windows <- seq(1, 91, by = window_size)
#     windowDist <- matrix(rep(-1, (length(windowStarts)-1)^2), nrow = length(windowStarts)-1)
#     windowMatchesList_Names <- vector("list", length(sam_mat2))
#     windowMatchesList_Nums <- vector("list", length(sam_mat2))
#     windowDistMatList <- vector("list", length(sam_mat2))
#     for(l in 1:length(sam_mat2)){
#         windowsMatchesNames <- vector("list", length(windowStarts)-1)
#         windowsMatchesNums <- vector("list", length(windowStarts)-1)
#
#         for(i in 1:(length(windowStarts)-1)){
#             windowsMatchesNames[[i]] <- dim_names[vapply(windowStarts[i]:windowEnds[i],
#                                                 function(x){
#                                                     which.max(sam_mat2[[l]][,x])
#                                                     }, numeric(1))]
#             windowsMatchesNums[[i]] <- vapply(windowStarts[i]:windowEnds[i],
#                                               function(x){
#                                                   max(sam_mat2[[l]][,x])
#                                                   }, numeric(1))
#
#         }
#         windowMatchesList_Names[[l]] <- windowsMatchesNames
#         windowMatchesList_Nums[[l]] <- windowsMatchesNums
#     }
#
#     newWindowMatchesList <- windowMatchesList_Nums
#
#     matchIdx <- lapply(1:length(windowMatchesList_Nums), function(x){
#         threshold <- tail(head(sort(as.vector(sam_mat2[[x]]),
#                                     decreasing = TRUE), 10),1)
#         lapply(windowMatchesList_Nums[[x]], function(y){which(y > threshold)})
#     })
#
#     ## regularize per factor with their individual thresholds
#     sam_mat_thresh <- sam_mat
#     for(i in 1:ncol(sam_mat)){
#         threshold <- tail(head(sort(sam_mat_thresh[,i], decreasing = TRUE), 10),1)
#         sam_mat_thresh[ which(sam_mat_thresh[,i] < threshold) , i] <- 0.0
#     }
#     heatmap3::heatmap3(sam_mat, Rowv = NA, Colv = NA, scale = "none",
#                        col = c("white", "yellow", "darkgreen", "darkblue", "orange"),
#                        revC = TRUE)
#     heatmap3::heatmap3(sam_mat_thresh, Rowv = NA, Colv = NA, scale = "none",
#                        col = c("white", "yellow", "darkgreen", "darkblue", "orange"),
#                        revC = TRUE)
#
#     sam_mat2_thresh <- lapply(sam_mat2,
#                               function(x){
#                                     threshold <- tail(head(sort(as.vector(x), decreasing = TRUE), 10),1)
#                                     x[x < threshold] <- 0.0
#                                     x
#                                 })
#     heatmap3::heatmap3(sam_mat2[[29]], Rowv = NA, Colv = NA, scale = "none",
#                        #col = c("white", "yellow", "darkgreen", "darkblue", "orange"),
#                        revC = TRUE)
#     heatmap3::heatmap3(sam_mat2_thresh[[29]], Rowv = NA, Colv = NA, scale = "none",
#                        #col = c("white", "yellow", "darkgreen", "darkblue", "orange"),
#                        revC = TRUE)
#     ##
#     sam_mat_thresh <- apply(sam_mat, MARGIN = 2,
#                             function(x){
#                                 threshold <- tail(head(sort(x, decreasing = TRUE), 10),1)
#                                 x[which(x < threshold)] <- 0.0
#                             })
#     ##
#
#     newWindowMatchesList_Nums <- lapply(1:length(windowMatchesList_Nums),
#                                         function(z){
#                                             lapply(1:length(windowMatchesList_Nums[[z]]),
#                                                    function(y){
#                                                        windowMatchesList_Nums[[z]][[y]][matchIdx[[z]][[y]]]})
#                                         })
#
#     newWindowMatchesList_Names <- lapply(1:length(windowMatchesList_Names),
#                                          function(z){
#                                              lapply(1:length(windowMatchesList_Names[[z]]),
#                                                     function(y){
#                                                         windowMatchesList_Names[[z]][[y]][matchIdx[[z]][[y]]]})
#                                          })
#
#     for(l in 1:length(windowDistMatList)){
#         a <- 58
#         b <- 64
#         for(i in 1:nrow(windowDist)){
#             for(j in 1:ncol(windowDist)){
#                 windowDist[i,j] <- length(intersect(newWindowMatchesList_Names[[a]][[i]],
#                                                     newWindowMatchesList_Names[[b]][[j]]))
#             }
#         }
#         # windowDistMatList[[l]] <- list()
#     }
#
#     a <- 58; b <- 64
#     for(i in 1:nrow(windowDist)){
#         for(j in 1:ncol(windowDist)){
#             windowDist[i,j] <- length(intersect(newWindowMatchesList_Names[[a]][[i]],
#                                                 newWindowMatchesList_Names[[b]][[j]]))
#         }
#     }
#     heatmap3::heatmap3(windowDist, Rowv = NA, Colv = NA, scale = "none",
#                        col = c("white", "yellow", "darkgreen", "darkblue", "orange"),
#                        revC = TRUE)
#     ##
#     # a <- 58; b <- 64
#     # for(i in 1:nrow(windowDist)){
#     #     for(j in 1:ncol(windowDist)){
#     #         windowDist[i,j] <- newWindowMatchesList_Nums[[a]][[i]],
#     #                                             newWindowMatchesList_Nums[[b]][[j]]
#     #     }
#     # }
#     # heatmap3::heatmap3(windowDist, Rowv = NA, Colv = NA, scale = "none",
#     #                    col = c("white", "yellow", "darkgreen", "darkblue", "orange"),
#     #                    revC = TRUE)
#
#
#     ##
# }


# .pos_agnostic_clust2 <- function(basisMat2,
#                                  positions,
#                                  normalize = TRUE,
#                                  windowSize = c(5,5)){
#     matAsVec <- apply(basisMat2, 2, function(y){y/max(y)})
#     print("done")
#     basisMat <- matrix(matAsVec,
#                        nrow = nrow(basisMat2), byrow = FALSE)
#     print(max(basisMat))
#
#     matchIdx <- lapply(1:length(windowMatchesList_Nums), function(x){
#         threshold <- tail(head(sort(as.vector(sam_mat2[[x]]),
#                                     decreasing = TRUE), 10),1)
#         lapply(windowMatchesList_Nums[[x]], function(y){which(y > threshold)})
#     })
#
# }


# .position_agnostic_clustering <- function(basisMat2,
#                                          positions,
#                                          topN = 10,
#                                          normalize = TRUE,
#                                          windowSize = c(5,5)){
#
#     # ## Considering using this if the scores are used in counting intersection
#     # startTime <- Sys.time()
#     #
#     ## normalize
#     matAsVec <- apply(basisMat2, 2, function(y){y/max(y)})
#     print("done")
#     ## put it back as a matrix
#     basisMat <- matrix(matAsVec, nrow = nrow(basisMat2), byrow = FALSE)
#     print(max(basisMat))
#
#
#     leftFlank <- windowSize[1]
#     rightFlank <- windowSize[2]
#     ## Starts and Ends of positions corresponding to the various dimers are
#     ##
#     borderEnds <- seq(length(positions), nrow(basisMat), by=length(positions))
#     borderStarts <- seq(1, nrow(basisMat), by=length(positions))
#
#     ##
#     ## impPos include the top-N positions/dimers -
#     ## where is top-10 threshold computed?
#     ## Values in this variable are numbers in the concatenated vector
#     ## (concatanetion of vectors denoting dimer presence, per dimer.
#     ## Thus, values here range from [1, 16*length(sequence)])
#     ##
#     impPos_dimerPosComb <- lapply(seq_len(ncol(basisMat)), function(x){
#         threshold <- tail(head(sort(as.vector(basisMat[,x]),
#                                     decreasing = TRUE), topN),1)
#         which(basisMat[,x] >= threshold)
#         })
#
#     ## Values in this variable are restricted to range [1, length(sequence)]
#     impPos_PosInSeq <- lapply(impPos_dimerPosComb,
#                               function(x){
#                                   tempAns <- x %% length(positions)
#                                   if(any(tempAns == 0)){
#                                       tempAns[which(tempAns == 0)] <- 1
#                                   }
#                                   return(tempAns)
#                                   })
#
#     impPos_PosInSeq_sorted <- lapply(impPos_PosInSeq, sort)
#     impPos_PosInSeq_sorted_idx <- lapply(impPos_PosInSeq, function(x){
#         sort(x, index.return = TRUE)$ix })
#
#     asMotifs <- lapply(impPos_PosInSeq_sorted, function(x) {
#         find_contig_motifs}
#         )
#
#
#     print(length(impPos_dimerPosComb))
#     print(impPos_dimerPosComb)
#     ##
#     leftFlankImpPos <- lapply(impPos, function(x){
#         x - leftFlank
#         })
#     rightFlankImpPos <- lapply(impPos, function(x){
#         x + rightFlank
#         })
#     ## dimerBin -- this gives the dimer ID (from 1 to 16) ofthe important positions
#     dimerBin <- lapply(impPos_dimerPosComb, function(x){ 1+floor(x/length(positions)) })
#     ## Same dimerBin range denotes the same dimer satifyign the threshold
#     ## For each factor, we should compare shifts for all those which have the same dimerBin
#     dimers <- get_dimers_from_alphabet(c("A", "C", "G", "T"))
#
#     dimerBin_sorted <- lapply(seq_along(dimerBin), function(x){
#         dimerBin[[x]][ impPos_PosInSeq_sorted_idx[[x]] ]
#     })
#
#     intersectCountP <- matrix(rep(0, ncol(basisMat)*ncol(basisMat)), nrow = ncol(basisMat))
#     # intersectCountN <- matrix(rep(0, ncol(basisMat)*ncol(basisMat)), nrow = ncol(basisMat))
#     ##
#     ## Motif-based comparison
#     for(i in 1:length(asMotifs)){
#         for(j in i:length(asMotifs)){
#             ## iterate over motifs of i
#             for(k in asMotifs[[i]]$starts){
#
#             }
#         }
#     }
#
#
#     # for(i in 1:length(impPos)){
#     #     a <- impPos[[i]]
#     #     print(a)
#     #     for(j in i:length(impPos)){
#     #         b <- impPos[[j]]
#     #         ## Make list of motifs -- either single or multiple dimers at
#     #         ## consecutive positions
#     #         # aList <- find_contig_motifs(a)
#     #         # bList <- find_contig_motifs(b)
#     #
#     #
#     #         ## Comparing all dimers in both factors
#     #         for(y in 1:length(a)){
#     #             for(z in 1:length(b)){
#     #                 if( dimerBin[[i]][y] == dimerBin[[j]][z] ){
#     #                     ## dimers match
#     #                     # if(leftFlankImpPos[[i]][y] < rightFlankImpPos[[j]][z] &&
#     #                     #    rightFlankImpPos[[i]][y] > leftFlankImpPos[[j]][z] ){
#     #                         ## shift only one of the pair
#     #                     if(leftFlankImpPos[[i]][y] <= impPos[[j]][z] &&
#     #                        impPos[[j]][z]<= rightFlankImpPos[[i]][y]){
#     #                         ## Overlap // good
#     #                         intersectCountP[i,j] <- intersectCountP[i,j] + 1
#     #                         #
#     #                         #*basisMat[impPos[[i]][y],i]*basisMat[impPos[[j]][z],j]
#     #                     }else{
#     #                         ## No Overlap // bad
#     #                         # intersectCountN[i,j] <-
#     #                         #     intersectCountN[i,j] -
#     #                         #     1*basisMat[impPos[[i]][y],i]*basisMat[impPos[[j]][z],j]
#     #                     }
#     #                 }else{
#     #                     ## dimers do not match // bad
#     #                     # intersectCountN[i,j] <-
#     #                     #     intersectCountN[i,j] -
#     #                     #     1*basisMat[impPos[[i]][y],i]*basisMat[impPos[[j]][z],j]
#     #                 }
#     #             }
#     #         }
#     #
#     #         # ## Only comparing selected dimers, i.e.,
#     #         # ## which dimers are present in both
#     #         # dimerMatch <- intersect(dimerBin[[i]], dimerBin[[j]])
#     #         # # message("Matching dimers: ", paste0(dimers[dimerMatch], sep=","))
#     #         # for(x in dimerMatch){
#     #         #     # print(paste0("Dimer: ",x ))
#     #         #     matchIdxA <- which(dimerBin[[i]] == x)
#     #         #     matchIdxB <- which(dimerBin[[j]] == x)
#     #         #     # message("Matches lengths: ", length(matchIdxA), "-",  length(matchIdxB))
#     #         #     ## check overlaps at these indices in the windows of important positions
#     #         #     for(y in matchIdxA){
#     #         #         for(z in matchIdxB){
#     #         #             ## check overlap, update count iff overlapping
#     #         #             if(leftFlankImpPos[[i]][y] < rightFlankImpPos[[j]][z] &&
#     #         #                rightFlankImpPos[[i]][y] > leftFlankImpPos[[j]][z] ){
#     #         #                 # message("Match/Overlap: ", paste0(c(i, j), sep="-" ) )
#     #         #                     intersectCountP[i,j] <-
#     #         #                         intersectCountP[i,j] +
#     #         #                         1*basisMat[impPos[[i]][y],i]*basisMat[impPos[[j]][z],j]
#     #         #             }else{
#     #         #                 # message("No Overlap: ", paste0(c(i, j), sep="-" ) )
#     #         #             }
#     #         #         }
#     #         #     }
#     #         #     ## normalize by:
#     #         #     ## how diverse in the matching dimer set (many distinct dimers are better)
#     #         #     ##
#     #         # }
#     #     }
#     # }
#
#     if(normalize){
#         matAsVec <- apply(intersectCountP, 2, function(y){y/max(y)})
#         print("normalize to max. 1")
#         intersectCountP <- matrix(matAsVec,
#                            nrow = nrow(intersectCountP), byrow = FALSE)
#     }
#     intersectCountP
#
# }


find_contig_motifs <- function(vector_positions){
    lagged_diffs <- diff(vector_positions)
    starts <- c()
    ends <- c()
    foundStart <- FALSE
    for( i in 1:length(lagged_diffs)){
        # print(i)
        if(!foundStart && lagged_diffs[i] == 1){
            starts <- append(starts, i)
            foundStart <- TRUE
        }
        else if(foundStart && lagged_diffs[i] == 1){
            ## do nothing
        }
        else if(foundStart && lagged_diffs[i] != 1){
            foundStart <- FALSE
            ends <- append(ends, i)
        }
        if(i == length(lagged_diffs) &&
           lagged_diffs[i] == 1 &&
           foundStart){
            ends <- append(ends, i+1)
        }
    }

    notOnes <- which(lagged_diffs != 1)
    toBeIncl <- setdiff(notOnes, ends)

    starts <- append(starts, toBeIncl)
    ends <- append(ends, toBeIncl)
    ## check starts and ends are of same length
    stopifnot(length(starts) == length(ends))
    ## check if last position is present in ends
    ## (it gets left out when not part of a contig)
    if(!any(ends == length(vector_positions))){
        starts <- append(starts, length(vector_positions))
        ends <- append(ends, length(vector_positions))
    }
    starts <- sort(starts)
    ends <- sort(ends)
    # print(paste(starts, ends, sep="-", collapse= ", "))
    motifs <- list(starts = starts,
                   ends = ends)
}




#' #' Keep this internal or external function?
#' #' @description We use hclust/hierarchical clustering for reordering archR
#' #' clusters
#' #'
#' #' @title Reorder archR clustering at the current last level
#' #'
#' #' @param archRresult The archRresult object
#' #'
#' #' @return Ordering returned from hclust
#' #'
#' #' @importFrom stats hclust dist
#' #' @export
#' reorder_archRresult <- function(archRresult, iteration = 3,
#'                                 clustMethod = "hc",
#'                                 linkage = "average",
#'                                 distMethod = "euclid",
#'                                 regularize = TRUE,
#'                                 topN = 10,
#'                                 returnOrder = FALSE,
#'                                 position_agnostic_dist = FALSE,
#'                                 config) {
#'     # Depends on archRresult object having a fixed set of names.
#'     # We need to .assert them
#'     # Finally, arrange clusters from processed outer chunks using hclust
#'     lastLevel <- iteration
#'
#'     basisMat <- archRresult$clustBasisVectors[[lastLevel]]$basisVector
#'
#'     if(regularize){
#'         basisMat2 <- basisMat
#'         for(i in 1:ncol(basisMat)){
#'             asVec <- as.vector(basisMat[,i])
#'             threshold <- tail(head(sort(asVec, decreasing = TRUE), topN),1)
#'             basisMat2[(basisMat[,i] < threshold), i] <- 0.0
#'         }
#'         basisMat <- basisMat2
#'     }
#'
#'
#'     if(position_agnostic_dist){
#'         #TODO
#'         factorsClustering <- .position_agnostic_clustering()
#'
#'     }else{
#'         setReturnOrder <- FALSE
#'         if(returnOrder){
#'             setReturnOrder <- TRUE
#'         }
#'         factorsClustering <- .handle_clustering_of_factors(basisMat,
#'                                       clustMethod = clustMethod,
#'                                       linkage = linkage,
#'                                       distMethod = distMethod,
#'                                       returnOrder = FALSE,
#'                                       flags = config$flags)
#'     }
#'     seqClusters <- get_seqs_clusters_in_a_list(archRresult$seqsClustLabels[[lastLevel]])
#'
#'     clusters <- .collate_clusters2(factorsClustering,
#'                        seqClusters)#, config$flags)
#'
#'     cluster_sol <- list(factorsClustering = factorsClustering,
#'                         clusters = clusters)
#'
#'     return(cluster_sol)
#'
#'
#'     # temp_hclust <-
#'     #     stats::hclust(stats::dist(
#'     #         t(archRresult$clustBasisVectors[[lastLevel]]$basisVectors)),
#'     #                          method = "ave")
#'     # new_order <- temp_hclust$order
#'     # # # Fetch original seqClustLabels as a list
#'     # origSeqsClustersAsList <- get_seqs_clusters_in_a_list(
#'     #                     archRresult$seqsClustLabelsList,
#'     #                     chooseLevel = lastLevel + 1)
#'     # if (length(new_order) != length(origSeqsClustersAsList)) {
#'     #     print("===== SAMARTH SAMARTH =====")
#'     #     print(archRresult$seqsClustLabels)
#'     #     print(origSeqsClustersAsList)
#'     #     print(length(origSeqsClustersAsList))
#'     #     print(new_order)
#'     #     print(length(new_order))
#'     #     stop("SAMARTH: Error")
#'     # }
#'     # ## arrange by the new ordering
#'     # newSeqsClusters <- lapply(new_order, function(x){
#'     #                                         origSeqsClustersAsList[[x]]
#'     #                                         }
#'     #                         )
#'     # ## the labels should be set as.character
#'     # newSeqsClustLabels <- unlist(lapply(
#'     #                                 seq_along(newSeqsClusters),
#'     #                                 function(x){
#'     #                                     as.character(
#'     #                                     rep(x, length(newSeqsClusters[[x]]))
#'     #                                     )
#'     #                                 }
#'     #                             )
#'     #                         )
#'     # ##
#'     # newClustBasisVectors <-
#'     #     archRresult$clustBasisVectors[[lastLevel]]$basisVectors[, new_order]
#'     # ##
#'     # new_field <- list(seqsClusters = newSeqsClusters,
#'     #                     seqsClustLabels = newSeqsClustLabels,
#'     #                     clustBasisVectors = newClustBasisVectors)
#'     # ##
#'     # archRresult$final <- new_field
#'     # ##
#'     # return(archRresult)
#' }
## =============================================================================

#' @title Retrieve sequence clusters as a list
#' @description Given the sequence cluster labels from a archR result object,
#' returns the clusters separated as a list.
#'
#' @param seqsClustLabels Sequences with cluster labels as in the archR result
#' object
#' @param chooseLevel Specify the level (archR iteration) at which sequence
#' clusters are to be reported. Default is 1.
#'
#' @return A list holding sequence clusters
#' @export
get_seqs_clusters_in_a_list <- function(seqsClustLabels, chooseLevel = 1){
    ## TODO: check that labels are not empty/NULL
    chosenLevelLabels <- seqsClustLabels
    clusterLevels <- levels(as.factor(seqsClustLabels))

    # chosenLevelLabels <- collect_cluster_labels(seqsClustLabels,
    #                                             chooseLevel = chooseLevel)
    # clusterLevels <- levels(as.factor(chosenLevelLabels))
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
