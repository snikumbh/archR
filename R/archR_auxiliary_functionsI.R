# @title Get cluster labels for chosen level/iteration
#
# @description Given a seqsClustLabels, collect cluster labels for sequences at
#  the chosen iteration/level.
#
# @param given_seqsClustLabels from archR result object
# @param chooseLevel choose a level/iteration. This value is number of
# iterations + 1
#
# @return A numeric vector of the same size as seqsClustLabels, with labels
# only up to the chosen iteration
#'
# collect_cluster_labels <- function(given_seqsClustLabels, chooseLevel = 1) {
#     ## Check if all_ok, all elements should have same length
#     .assert_archR_seqsClustLabels_at_end(given_seqsClustLabels)
#     splitChar <- "-"
#     elements_length <- unique(unlist(lapply(strsplit(given_seqsClustLabels,
#                                                     split = splitChar),
#                                                     length)))
#     if (chooseLevel > elements_length) {
#         stop("Choose_levels (", chooseLevel,
#         ") greater than levels present in cluster_labels(",
#         elements_length, ").")
#     } else {
#         selectedLabels <- unlist(lapply(strsplit(given_seqsClustLabels,
#                                                 split = splitChar),
#                                         function(x) {
#                                             paste0(x[seq_len(chooseLevel)],
#                                                 collapse = splitChar)
#                                         }))
#     }
#     .assert_archR_seqsClustLabels_at_end(selectedLabels)
#     return(selectedLabels)
# }
## =============================================================================

# @title Assign samples to clusters
#
# @param samplesMatrix
#
# @return A list that can be assigned as an element in globClustAssignments
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



# .get_consensusClustMemberships <- function(allRunMemberships, nClusters) {
# 
#     consensusClustMemberships <-
#         apply(allRunMemberships, 2, function(yVec){
#                         which.max(
#                             vapply(seq_len(nClusters),
#                                    function(x){sum(yVec == x)},
#                             numeric(1))
#                         )}
#         )
# 
#     return(consensusClustMemberships)
# }


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
# .map_clusters_to_factors <- function(samplesMatrix, clustOrderIdx, iChunksColl,
#                                     iChunkIdx, flags) {
#     ## clustOrderIdx IS clustering_sol_kmeans$reordering_idx, the reorderingIdx
#     ## is a nested list, i.e. a list holding list of seqs_ids belonging to one
#     ## cluster together. Hence, length(this variable) gives #clusters.
#     ## samplesMatrix IS samplesMatrix.
#     ##
#     ## Make checks
#     .assert_archR_samplesMatrix(samplesMatrix)
#     ##
#     if (!is.list(clustOrderIdx)) {
#         print(clustOrderIdx)
#         stop("Mapping clusters to factors, cluster orders is not a list")
#     } else if (!is.vector(clustOrderIdx[[1]])) {
#             stop("Mapping clusters to factors, cluster orders is not a vector")
#     }
#     ## check iChunksColl and iChunkIdx
#     if (!is.null(iChunksColl) && length(iChunksColl) < 1) {
#         stop("Mapping clusters to factors, expecting inner chunks as a
#                 list of length > 0")
#     }
#     if (!is.null(iChunkIdx) && !is.numeric(iChunkIdx) && iChunkIdx > 0) {
#         stop("Mapping clusters to factors, expecting inner chunk index as a
#                 numeric > 0")
#     }
#     ##
#     .assert_archR_flags(flags)
#     ##
#     rightClusterOrders <- vector("list", length(clustOrderIdx))
#     if (length(clustOrderIdx) == 1) {
#         ## Special case
#         rightClusterOrders[[1]] <- iChunksColl[[iChunkIdx]][clustOrderIdx[[1]]]
#         return(rightClusterOrders)
#         ##
#     } else if (length(clustOrderIdx) > 1) {
#         ##
#         message("Cluster to Factor:")
#         for (cluster_idx in seq_along(clustOrderIdx)) {
#             ##
#             # print(apply(as.matrix(samplesMatrix[, clustOrderIdx[[cluster_idx]]]),
#             #             1, quantile))
#             # print(as.matrix(samplesMatrix[, clustOrderIdx[[cluster_idx]]]))
#             rowwiseMeans <- rowMeans(
#                 as.matrix(samplesMatrix[, clustOrderIdx[[cluster_idx]]])
#             )
#             print(paste0("Cluster ID:", cluster_idx))
#             print(paste0("Rowwise means", rowwiseMeans))
#             relevant_factor <- which.max(rowwiseMeans)
#             relevant_descending_order <-
#                 sort(samplesMatrix[relevant_factor,clustOrderIdx[[cluster_idx]]],
#                      index.return = TRUE,
#                      decreasing = TRUE)$ix
#             # print(relevant_factor)
#             # #stop("SAMARTH")
#             ##
#             # relevant_factor <-
#             #     which.max(rowMeans(
#             #     as.matrix(samplesMatrix[, clustOrderIdx[[cluster_idx]]])
#             #     ))
#             if (length(relevant_factor) > 1) {
#                 if (flags$debugFlag) {
#                     message("Multiple basis vectors attained max.")
#                     message(relevant_factor)
#                 }
#             } else {
#                 message(relevant_factor)
#                 if (flags$debugFlag) {
#                     message("Factor-mapping OK")
#                     message(relevant_factor)
#                 }
# 
#                 seqs_vector_to_assign <-
#                     iChunksColl[[iChunkIdx]][clustOrderIdx[[cluster_idx]]]
#                 ## order the sequences by the scores (highest first, weakest last)
#                 rightClusterOrders[[relevant_factor]] <-
#                     seqs_vector_to_assign[relevant_descending_order]
#             }
#         }
#         ## Check if any factor got left out?
#         ## i.e. no seqs assigned to its cluster
#         if (any(lapply(rightClusterOrders, length) == 0)) {
#             print(clustOrderIdx)
#             thisGotLeftOut <- which(lapply(rightClusterOrders, length) == 0)
#             warning(c("Factor(s) got no sequences assigned: ",
#                     paste0(thisGotLeftOut, sep="-")), immediate. = TRUE)
#             for (id in thisGotLeftOut) {
#                 rightClusterOrders[[id]] <- 0
#             }
#         }
#         return(rightClusterOrders)
#     }
# }
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
    ## Important: The collatedClustAssignments variable is a list where each 
    ## element holds the indices of sequences falling in the same cluster, like so
    ## For a set of 200 sequences, collatedClustAssignments is
    ## [[1]]
    ## [1]   1   2   6  11  20  23  32  46  50  52  63  68  71  72  75  80  82  83  85  86  87  88  92  99 102 104 105
    ## [28] 112 114 116 123 131 137 138 143 145 146 147 148 151 158 165 167 168 171 176 179 182 190 197
    ##
    ## [[2]]
    ## [1]   3   5   7   8  10  12  13  16  17  18  19  22  25  30  31  36  37  38  39  41  43  45  47  49  51  53  54
    ## [28]  55  56  58  59  60  61  64  67  69  70  73  77  84  89  91  94  95  98 103 106 107 109 115 117 118 119 120
    ## [55] 121 122 125 126 128 130 132 133 134 135 136 140 141 142 144 150 152 153 154 156 159 160 161 162 164 166 173
    ## [82] 175 177 180 183 184 185 186 188 189 191 193 194 195 196 198 199 200
    ##
    ## [[3]]
    ## [1]   4   9  14  15  21  24  26  27  28  29  33  34  35  40  42  44  48  57  62  65  66  74  76  78  79  81  90
    ## [28]  93  96  97 100 101 108 110 111 113 124 127 129 139 149 155 157 163 169 170 172 174 178 181 187 192
    ## 
    ## This is different than the globClustAssignment variable which will hold
    ## the indices of the factors that are combined into one cluster. For instance,
    ## when there are 25 factors/clusters which have combined to give 5 clusters,
    ## globClustAssignments will hold something like this:
    ## [[1]]
    ## [1] "3"  "13" "11" "16" "6" 
    ##
    ## [[2]]
    ## [1] "15" "1"  "8"  "25" "18"
    ##
    ## [[3]]
    ## [1] "24" "5"  "4"  "10" "21"
    ## 
    ## [[4]]
    ## [1] "14" "12" "17" "2"  "22"
    ## 
    ## [[5]]
    ## [1] "19" "9"  "20" "7"  "23"
    ## 
    ## The collatedClustAssignments variable and the globClustAssignments variable were confused,
    ## and a corresponding test was wrongly written marking the previous update_cluster_labels code block 
    ## as erroneous. However, that is not correct, as has been noted above.
    ## For reference: see commit 4705553
    ## 
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
    newSeqsClustLabels <- oldSeqsClustLabels
    for (i in seq_len(nClusters)) {
        needUpdateIdx <- collatedClustAssignments[[i]]
        newSeqsClustLabels[needUpdateIdx] <-
            vapply(newSeqsClustLabels[needUpdateIdx],
                    function(x) {
                        paste0(candidateClustLabels[i])
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


#' @title Retrieve sequence clusters as a list
#' @description Given the sequence cluster labels from a archR result object,
#' returns the clusters separated as a list.
#'
#' @param seqsClustLabels Sequences with cluster labels as in the archR result
#' object.
#'
#' @return A list holding sequence clusters.
#' @export
get_seqs_clusters_in_a_list <- function(seqsClustLabels){
    ## check that labels are not empty/NULL
    .assert_archR_seqsClustLabels(seqsClustLabels)
    clusterLevels <- levels(as.factor(seqsClustLabels))

    seqs_clusters_as_a_list <- lapply(seq_along(clusterLevels),
                                        function(x){
                                            which(seqsClustLabels ==
                                                    clusterLevels[x])
                                    })

    return(seqs_clusters_as_a_list)
}
## =============================================================================
