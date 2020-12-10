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
    ## When chunkLength (i.e., total sequences) < reqdChunkSize, the else 
    ## condition, return the totalSet as the only chunk.
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
                                        returnOrder = TRUE,
                                        useCutree = TRUE,
                                        parentChunks = NULL) {
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
        if(flags$debugFlag) message("Hierarchical clustering")
        ############### REORDER CLUSTERS FROM OUTER CHUNK ######################
        if(flags$debugFlag) message("RE-ORDERING CLUSTERS")

        if(flags$debugFlag) {
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
            if(useCutree){
                clustList <- .get_clusters_from_hc_using_cutree(
                                    hcObj = temp_hclust,
                                    distMat = as_dist_mat,
                                    hStep = 0.05, 
                                    parentChunks = parentChunks,
                                    verbose = flags$debugFlag)
            }else{
                clustList <- .get_clusters_from_hc(hcObj = temp_hclust,
                                           distMat = globFactorsDistMat,
                                           verbose = flags$debugFlag)
            }
            if(length(clustList) == 1){
                clustList <- lapply(seq_len(max(clustList[[1]])), 
                                    function(x){x})   
            }
            if(flags$debugFlag) {
                message("=== DONE ===")
                print(clustList)
            }
            return(clustList)
        }
    }
}
## =============================================================================

################################# IMPROVED VERSION #############################

#' @title Get clusters from hclust object (improved version of previous ad-hoc
#' approach -- this allows singletons)
#' 
#' @param hcObj The hierachical clustering, hclust object
#' @param linkage String. Specify the linkage used -- works best for ward.D.
#' @param bound_scale Specify a threshold percentage. This will be used to 
#' compute the extent of the distance allowed between any two elements so that
#' they can be combined. See more explanation in the function code
#' @param pIter Numebr of iterations of joining procedures to be performed (See
#' while loop in the function code).
#' 
#' @importFrom utils tail
#' @importFrom stats median
.get_clusters_from_hc2 <- function(hcObj, linkage = "complete", 
                                   bound_scale = 0.75, pIter = 4,
                                   distMat, verbose = FALSE){
    ## important notes:
    ## -- param bound_scale: changing the value of bound_scale affects formation
    ##    of singleton clusters
    ## -- param linkage: currently only works for/accepts 'complete' linkage. 
    ##    Average linkage leaves many unjoined in remainingIdx
    
    
    #####
    ## Algo to cut dendrogram and form clusters
    ## corner case: 1. what if start nodes are not available?
    ##
    ## Improvements over the previous version:
    ## 1. Form start nodes using the merge information. These start nodes are 
    ## pairs of indices that were first merged when forming the dendrogram.
    ## 2. Join additional nodes 
    ## 
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
    
    ## Assign start nodes
    ## Use merge information from dendrogram object to form start nodes of 
    ## paired indices
    startNodeList <- list()
    ## collect rows from hcObj$merge where both indices col1 and col2 are < 0
    relRowIdx_from_merge <- which(hcObj$merge[,1] < 0 & hcObj$merge[,2] < 0)
    ## Are all pair-merges equivalent? Some pair-merges could happen at a much larger 
    ## distance than others. Worth detecting and not adding them?
    # pairMergesDist <- unlist(lapply(relRowIdx_from_merge, function(x){
    #     c1 <- -1*hcObj$merge[x,1]
    #     c2 <- -1*hcObj$merge[x,2]
    #     distMat[c1,c2]
    # }))
    # median(pairMergesDist) - 3*stats::mad(pairMergesDist)
    # median(pairMergesDist) + 3*stats::mad(pairMergesDist)
    ##
    startNodeList <- lapply(relRowIdx_from_merge, function(x){
        -1*hcObj$merge[x,]
    })
    ## ^ start nodes formed
    
    ## Grow start nodes by joining neighboring indices from hcObj$order
    
    nodeList <- startNodeList
    remainingIdx <- setdiff(hcObj$order, unlist(nodeList))
    remIter <- pIter ## 
    while(length(remainingIdx) < 1 || remIter > 0){
        ## Looking to join left neighbors
        ## -- Check which node in nodeList has a candidate left neighbor to be 
        ##    considered for stringing/joining
        ## -- Decide if this candidate left neighbor will be joined based on a 
        ##    criterion: (a) it is not already included in the nodeList
        ##               (b) dist between the neighbor and the elements of the 
        ##                   node is too high (increases the within-cluster 
        ##                   some measure)
        ## --            (c) Is it close enough this cluster (of start nodes) to 
        ##                   justify joining
        nodeList <- lapply(nodeList, function(x){
            if(verbose) print(c("Next:", paste0(x, collapse="-")))
            leftNeighborOfX <- NA
            leftNeighborOfXisFirst <- FALSE
            positionOfXinOrder <- which(x[1] == hcObj$order)
            # if(verbose) print("hcObj$order for ref:");print(hcObj$order)
            # if(verbose) print(paste("position in Order:", positionOfXinOrder))
            if(positionOfXinOrder != 1){ 
                ## not the first, so safely subtract 1 from index
                leftNeighborOfX <- hcObj$order[positionOfXinOrder-1]
                if(verbose) print(paste("Left neighbor:", leftNeighborOfX))
            }else{
                ## ## Else, it is the first node, without a 
                ## leftNeighborCandidate. Let it be NA
            }
            
            if(!is.na(leftNeighborOfX)){
                ## Already note if the leftNeighborOfX is the first node
                if(positionOfXinOrder == 2) leftNeighborOfXisFirst <- TRUE
                
                ## Decide if we can join leftNeighborOfX to node with X?
                decisionToJoinTrue <- TRUE
                ## There are conditions in which the decision to join should FALSE
                ##
                ## Case 1: leftNeighbor already in another node among nodeList
                if(any(leftNeighborOfX == unlist(nodeList))){
                    ## already in another node; so skip joining
                    decisionToJoinTrue <- FALSE
                    if(verbose) print("== Not left joined, case 1 ==")
                }
                ##
                ## Case 2: dist b/w leftNeighbor and x[1] is greater than 
                ##         leftNeighbor and leftNeighbor of leftNeighbor
                ##      -- This will hold for complete linkage, but not for 
                ##         average linkage.
                if(decisionToJoinTrue && !leftNeighborOfXisFirst){
                    ## can use dist to check if it cannot be joined
                    positionOfLeftOfLeft <- positionOfXinOrder - 2
                    leftOfLeftNeighbor <- hcObj$order[positionOfLeftOfLeft]
                    ##
                    # if(linkage != "average"){
                        if(distMat[x[1], leftNeighborOfX] > 
                                distMat[leftNeighborOfX, leftOfLeftNeighbor]){
                            decisionToJoinTrue <- FALSE
                            if(verbose) print(paste("== Not left joined, case 2", 
                                                      " (may join w/ other) ==", 
                                                    collapse = ""))
                        }
                    # }else{
                    #     if(distMat[x[1], leftNeighborOfX] > 
                    #        distMat[leftNeighborOfX, leftOfLeftNeighbor]){
                    #         decisionToJoinTrue <- FALSE
                    #         if(verbose) print(paste("== Not left joined, case 2", 
                    #                                 " (may join w/ other) ==", 
                    #                                 collapse = ""))
                    #     }else{
                    #         ## Add as a different node?
                    #         if(verbose) print("avg SAMARTH separate CLUSTER")
                    #         return(list(x, leftNeighborOfX))
                    #     }
                    # }
                }
                ##
                ## Case 3: LeftNeighbor should be kept as a separate node or 
                ##         effectively, a separate cluster
                ##
                idx_pairs <- t(utils::combn(x, 2))
                distInSet <- unlist(lapply(seq_len(nrow(idx_pairs)), function(y){
                    distMat[idx_pairs[y,1],idx_pairs[y,2]]
                })) 
                A <- stats::median(distInSet)
                dist_bound <- A + bound_scale*A
                ## If previous cases still kept decision to join as TRUE, check
                ## the distance criterion, and finally decide.
                ## If the previous cases already decided as FALSE, don't bother
                if(decisionToJoinTrue && 
                   distMat[leftNeighborOfX, x[1]] > dist_bound){
                    decisionToJoinTrue <- FALSE
                    if(verbose) print("== Not left joined, case 3 ==")
                    ## Add a separate node/cluster
                    if(verbose) print("SAMARTH separate CLUSTER")
                    return(list(x, leftNeighborOfX))
                }
                ## Exclude all above cases, then join
                if(decisionToJoinTrue){
                    ## If decision TRUE, join from the left (use after = 0)
                    x <- append(x, values = leftNeighborOfX, after = 0)
                    if(verbose){ print("joined x:"); print(x)}
                }
            }else{if(verbose)print("Already first, no left neighbor")}
            x
        })
        # ## unravel the nested list element
        nodeList <- unfurl_nodeList(nodeList)
        
        if(verbose){ print("Updated nodeList");print(paste(nodeList))}
        
        ## Update remainingIdx
        remainingIdx <- setdiff(hcObj$order, unlist(nodeList))
        if(length(remainingIdx) > 0){
            if(verbose){ print("What remains:");print(remainingIdx)}
        }else{
            if(verbose){ print("breaking loop:");print(remainingIdx)}
            break
        }
        if(verbose) message("l(remaining): ", length(remainingIdx))
        
        # if(linkage != "average"){
            if(verbose) message("=== Right joining ===")
            ## Join right side nodes
            ## Looking to join right neighbors
            ## -- Check which node in nodeList has a candidate right neighbor to be 
            ##    considered for stringing/joining
            ## -- Decide if this candidate right neighbor will be joined based on a 
            ##    criterion: (a) it is not already included in the nodeList
            ##               (b) dist between the neighbor and the elements of the 
            ##                   node is too high (increases some within-cluster 
            ##                   measure)
            ## --            (c) Is it close enough this cluster (of start nodes) to 
            ##                   justify joining 
            nodeList <- lapply(nodeList, function(x){
                if(verbose) print(c("Next:", paste0(x, collapse="-")))
                rightNeighborOfX <- NA
                rightNeighborOfXisLast <- FALSE
                totalLength <- length(hcObj$order)
                rightOfX <- utils::tail(x, 1)
                positionOfXinOrder <- which(rightOfX == hcObj$order)
                # if(verbose) print("hcObj$order for ref:");print(hcObj$order)
                # if(verbose) print(paste("position in Order:", positionOfXinOrder))
                if(positionOfXinOrder != totalLength){
                    ## not the last, so safely add 1 to index
                    rightNeighborOfX <- hcObj$order[positionOfXinOrder+1]
                    if(verbose) print(paste("Right neighbor:", rightNeighborOfX))
                }else{
                    ## Else, it is the last in order, without a 
                    ## rightNeighborCandidate. Let it be NA
                }
                
                if(!is.na(rightNeighborOfX)){
                    ## Already note if the rightNeighborOfX is the last node
                    if(positionOfXinOrder == totalLength-1) rightNeighborOfXisLast <- TRUE
                    
                    ## Decide if we can join rightNeighborOfX to node with X?
                    decisionToJoinTrue <- TRUE
                    ## There are conditions in which the decision to join should be FALSE
                    ##
                    ## Case 1: rightNeighbor already in another node among nodeList
                    if(any(rightNeighborOfX == unlist(nodeList))){
                        ## already in another node; so skip joining
                        decisionToJoinTrue <- FALSE
                        if(verbose) print("== Not left joined, case 1 ==")
                    }
                    ##
                    ## Case 2: dist b/w rightNeighbor and x[1] is greater than 
                    ##         rightNeighbor of rightNeighbor
                    if(decisionToJoinTrue && !rightNeighborOfXisLast){
                        ## can use dist to check if it cannot be joined
                        positionOfRightOfRight <- positionOfXinOrder + 2
                        rightOfRightNeighbor <- hcObj$order[positionOfRightOfRight]
                        ##
                        if(distMat[rightOfX, rightNeighborOfX] > 
                           distMat[rightNeighborOfX, rightOfRightNeighbor]){
                            decisionToJoinTrue <- FALSE
                            if(verbose) print(paste("== Not left joined, case 2", 
                            " (possibility to join with other ==", collapse = ""))
                        }
                    }
                    ##
                    ## Case 3: rightNeighbor should be kept as a separate node or 
                    ##         effectively, a separate cluster
                    ##
                    idx_pairs <- t(utils::combn(x, 2))
                    distInSet <- unlist(lapply(seq_len(nrow(idx_pairs)), function(y){
                        distMat[idx_pairs[y,1],idx_pairs[y,2]]
                    }))
                    dist_bound <- median(distInSet) + bound_scale*median(distInSet)
                    ## If previous cases still kept decision to join as TRUE, check
                    ## the distance criterion, and finally decide.
                    ## If the previous cases already decided as FALSE, don't bother
                    if(decisionToJoinTrue && 
                       distMat[rightNeighborOfX, rightOfX] > dist_bound){
                        decisionToJoinTrue <- FALSE
                        if(verbose) print("== Not left joined, case 3 ==")
                        ## Add a separate node/cluster
                        if(verbose) print("SAMARTH separate CLUSTER")
                        return(list(x, rightNeighborOfX))
                    }
                    ## Exclude all above cases, then join
                    if(decisionToJoinTrue){
                        ## If decision TRUE, join from the right (use after = length(x))
                        x <- append(x, values = rightNeighborOfX, after = length(x))
                        if(verbose){ print("joined x:"); print(x)}
                    }
                }else{if(verbose)print("Already last, no right neighbor")}
                ## return
                x
            })
            ########
            # ## unravel the nested list element
            nodeList <- unfurl_nodeList(nodeList)
            
            ## If all resolved already, break loop
            if(length(remainingIdx) == 0) break
            
            ########
            if(verbose){print("Updated nodeList");print(paste(nodeList))}
            
            ## Update remainingIdx
            remainingIdx <- setdiff(hcObj$order, unlist(nodeList))
            if(length(remainingIdx) > 0){
                if(verbose){ print("What remains:");print(remainingIdx)}
            }else{
                if(verbose){ print("breaking loop:");print(remainingIdx)}
                break
            }
        # }
        remIter <- remIter - 1
        if(verbose) message("remIter: ", remIter)
        if(verbose) message("l(remaining): ", length(remainingIdx))
    }
    
    ## sort the clusters based on hcObj$order
    starts <- lapply(nodeList, function(x){utils::head(x, 1)})
    matchOrder <- vapply(starts, function(x){which(x == hcObj$order)}, numeric(1))
    newIdx <- sort(matchOrder, index.return = TRUE)$ix
    nodeListUpd <- lapply(newIdx, function(x){nodeList[[x]]})
    nodeList <- nodeListUpd
    
    if(verbose) print(paste(nodeList))
    
    #####
    return(nodeList)
}

unfurl_nodeList <- function(nodeList){
    
    element_lengths <- unlist(lapply(nodeList, function(elem){
        if(class(elem) == "list") length(elem)
        else 1
    }))
    
    if(any(element_lengths != 1)){
        new_list <- vector("list", sum(element_lengths))
        iter1 <- 0
        iter2 <- 1
        while(iter2 <= length(nodeList)){
            if(class(nodeList[[iter2]]) != "list"){
                iter1 <- iter1 + 1
                new_list[[iter1]] <- nodeList[[iter2]]
            }else{
                for(i in seq_along(nodeList[[iter2]])){
                    iter1 <- iter1 + 1
                    new_list[[iter1]] <- nodeList[[iter2]][[i]]
                }
            }
            iter2 <- iter2 + 1
        }
        new_list
    }else{
        nodeList
    }
}
############################## IMPROVED VERSION ################################


############################## CUTREE VERSION ##################################
## This works best/is best referred with ward.D linkage
#' @title Get clusters out of hierarchical clustering object using cutree. 
#' Used internally by archR.
#' 
#' @description Clusters from a hierarchical clustering object are obtained 
#' by using cutree at different heights of the tree.
#' 
#' @param hcObj The hierarchical clustering object as returned by 
#' \code{\link{stats::hclust}}.
#' @param distMat The distance matrix that was used by hclust, a 
#' \code{\link{stats::dist}} object.
#' @param hStep Numeric. The step size used to increment height values for 
#' cutree. Default value is 0.05.
#' @param parentChunks List. Specify the factor numbers in the previous 
#' iteration of archR that factors in the current iteration resulted from. 
#' Default value is NULL.
#' @param verbose Logical. Specify TRUE for verbose output.
#' 
#' @importFrom stats cutree 
#' @importFrom cluster silhouette
#'
.get_clusters_from_hc_using_cutree <- function(hcObj, distMat, hStep = 0.05,
                                               parentChunks = NULL,
                                               verbose = FALSE){
    
    cut_heights <- seq(min(hcObj$height), max(hcObj$height), by = hStep)
    
    ## Important: Address the question of whether any merging is required?
    ## -- We could use information from archR iteration
    ##    -- Obtain the cutree-using-h clustering here. 
    ##    -- Is it combining subsequent factors, like 1-2-3, 3-4 etc.? 
    ##          a. This could mean that it is undoing the clusters identified by 
    ##          archR in a previous iteration.
    ##          b. But, if 1 came from chunk 1 and 2-3 came from chunk 2, and if
    ##          the clustering results in combining 1-2, then this is a plausible 
    ##          merge and should be left alone.
    ##          c. So, we use a parentChunks variable that notes the parent 
    ##          chunk for each ofthe current factors. This info can be used to 
    ##          make this decision unambiguously.
    
    if(verbose) startTm <- Sys.time()
    sils_cut_by_h <- unlist(lapply(cut_heights, 
                       function(x){
                           foo_try <- cutree(hcObj, h = x)
                           names(foo_try) <- NULL
                           if(length(unique(foo_try)) > 1){
                               sils <- cluster::silhouette(foo_try, 
                                                    dist = distMat)
                               base::mean(sils[, "sil_width"])
                           }else return(-100)
                       }))
    if(verbose) Sys.time() - startTm
    
    
    if(verbose) plot(sils_cut_by_h)
    
    ## multiple matches, first match index is returned with which.max
    cheight_idx <- which.max(sils_cut_by_h)
    
    
    if(verbose) message("Max. silhouette score at index = ", cheight_idx,
            ", height:", cut_heights[cheight_idx])
    
    cut_result <- stats::cutree(hcObj, h = cut_heights[cheight_idx])
    names(cut_result) <- NULL
    clust_list <- lapply(1:length(unique(cut_result)), 
                                     function(x) which(cut_result == x))
    
    if(verbose) message("#Clusters: ", length(clust_list))
    if(verbose) print(paste(clust_list))
    
    if(!is.null(parentChunks)){
        if(verbose) message("Checking parent chunks...")
        childrenPerParent <- lapply(unique(parentChunks), function(x){
                which(parentChunks == x)
            })
        updated_clust_list <- lapply(seq_along(clust_list), function(x){
                        parents <- unique(parentChunks[clust_list[[x]]])
                        thisX <- clust_list[[x]]
                        if(verbose) message("parents");print(parents)
                        if(verbose) message("thisX");print(thisX)
                        if(length(parents) == 1 && length(thisX) > 1){
                            ## if cluster elements come from same parent 
                            ## chunk of previous iteration
                            ## + now check if all elements in this cluster
                            ## are the only children of that parent chunk
                            ## (in other words nothing got separated and 
                            ## combined with other factors)
                            childrenThisParent <- 
                                childrenPerParent[[as.integer(parents)]]
                            if(verbose){
                                message("length is One so, childrenThisParent:")
                                print(childrenThisParent)
                                message(class(as.numeric(childrenThisParent)))
                                message(class(thisX))
                                message(class(as.numeric(thisX)))
                            }
                            ##
                            if(identical(as.numeric(childrenThisParent), 
                                         as.numeric(thisX))){
                                if(verbose) message("Are identical")
                                ## Split the merge back into separate clusters
                                return(lapply(childrenThisParent, function(x){x}))
                            }else{
                                if(verbose) print("samarth, returning thisX")
                                return(thisX)
                            }
                        }
                        else{
                            ## cluster elements come from different
                            ## parent chunks (of previous iteration)
                            return(thisX)
                        }
                    })
        clust_list <- unfurl_nodeList(updated_clust_list)
        if(verbose) message("Updated clusters")
    }
    
    if(verbose) message("#Clusters: ", length(clust_list))
    if(verbose) paste(clust_list)
    return(clust_list)
}

############################## CUTREE VERSION ##################################



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

    # if(verbose){
    #     message("Final Length: ", sum(unlist(lapply(nodeList, length))))
    # }
    stopifnot(sum(unlist(lapply(nodeList, length))) == length(newOrder))
    #####
    return(nodeList)
}


#' @title Reorder archR raw clustering at the chosen iteration
#' @description We use hierarchical clustering for reordering archR raw
#' clusters
#' @param archRresult The archRresult object.
#' @param iteration Specify clusters at which iteration of archR are to be reordered.
#' @param clustMethod Specify 'hc' for hierarchical clustering.
#' @param linkage One of linkage values as specified for hierarchical clustering.
#' @param distMethod Distance measure to be used with hierarchical clustering.
#' @param regularize Specify TRUE if regularization is to be performed before 
#' comparison. See argument 'topN'.
#' @param topN Keep only the topN positions for comparing basis vectors.
#' @param returnOrder Specify TRUE if only the computed order for hierarchical 
#' clustering is to be returned
#' @param position_agnostic_dist If position agnostic distance measure is to be 
#' used. For future implementation.
#' @param decisionToReorder Logical. Specify TRUE if reordering using 
#' hierachical agglomerative clustering is to be performed, otherwise FALSE. 
#' 
#' @param config Pass the configuration of the archR result object.
#'
#' @return Returns ordering returned from hclust or the re-ordered clusters. 
#' When `decisionToReorder' is FALSE, it returns the already existing basis 
#' vectors, each as singleton clusters. The sequence cluster labels and sequence
#'  clusters are handled accordingly.
#'
#' @importFrom stats hclust dist
#' @importFrom utils tail
#' @export
reorder_archRresult <- function(archRresult, iteration = 3,
                                clustMethod = "hc",
                                linkage = "average",
                                distMethod = "euclid",
                                regularize = TRUE,
                                topN = 10,
                                returnOrder = FALSE,
                                position_agnostic_dist = FALSE,
                                decisionToReorder = TRUE,
                                config) {
    # Depends on archRresult object having a fixed set of names.
    # We need to .assert them
    # Finally, arrange clusters from processed outer chunks using hclust
    
    basisMat <- archRresult$clustBasisVectors[[iteration]]$basisVector
    
    if(regularize){
        basisMat2 <- basisMat
        for(i in 1:ncol(basisMat)){
            asVec <- as.vector(basisMat[,i])
            threshold <- utils::tail(utils::head(sort(asVec, decreasing = TRUE), topN),1)
            basisMat2[(basisMat[,i] < threshold), i] <- 0.0
        }
        basisMat <- basisMat2
    }
    
    ## Prepare info on parent chunks for each cluster/factor at given iteration
    ## We get this info from the seqsClustLabels
    parentChunks <- NULL
    if(iteration > 1){
        clustsThisIter <- sort(unique(archRresult$seqsClustLabels[[iteration]]))
        parentChunks <- unlist(lapply(clustsThisIter, function(x){
            relevantSeqsIds <- which(archRresult$seqsClustLabels[[iteration]] 
                                     == x)
            unique(archRresult$seqsClustLabels[[iteration - 1]][relevantSeqsIds])
        }))
    }
    if(decisionToReorder){
        if(position_agnostic_dist){
            # TODO
            
        }else{
            setReturnOrder <- FALSE
            if(returnOrder){
                setReturnOrder <- TRUE
            }
            factorsClustering <- .handle_clustering_of_factors(basisMat,
                                               clustMethod = clustMethod,
                                               linkage = linkage,
                                               distMethod = distMethod,
                                               returnOrder = returnOrder,
                                               flags = config$flags,
                                               parentChunks = parentChunks)
            if(config$flags$debugFlag){
                message("Performed factor clustering, returning object:")
                print(factorsClustering)
            }
        }
        
    }else{
        ## Do not reorder, but prepare the return object
        nFactors <- archRresult$clustBasisVectors[[iteration]]$nBasisVectors
        factorsClustering <- vector("list", nFactors)
        factorsClustering <- lapply(1:nFactors, function(x){x})
        if(config$flags$debugFlag){
            message("No factor clustering, returning object:")
            print(factorsClustering)
        }
    }
    seqClusters <- get_seqs_clusters_in_a_list(archRresult$seqsClustLabels[[iteration]])
    
    clusters <- .collate_clusters2(factorsClustering,
                                   seqClusters)
    clustLabels <- .update_cluster_labels(oldSeqsClustLabels = 
                                              archRresult$seqsClustLabels[[iteration]], 
                                          collatedClustAssignments = clusters,
                                          flags = config$flags)
    
    cluster_sol <- list(basisVectorsClust = factorsClustering,
                        clusters = clusters,
                        seqsClustLabels = clustLabels)
    
    return(cluster_sol)

}
## =============================================================================


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
