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
# 

## Getter functions ============================================================

## get functions for components of archR's result object
## These keeps such accesses in other functions agnostic to how they are 
## maintained in the archR result object
get_seqClLab <- function(res, iter){
    
    if(is.null(iter)){
        return(res$clustSol$seqsClustLabels)
    }else{
        return(res$seqsClustLabels[[iter]])
    }
}

get_clBasVec_k <- function(res, iter){
    
    return(res$clustBasisVectors[[iter]]$nBasisVectors)
}

get_clBasVec_m <- function(res, iter){
    
    return(res$clustBasisVectors[[iter]]$basisVectors)
}

get_clBasVec <- function(res, iter){
    
    return(res$clustBasisVectors[[iter]])
}

#' @title seqs_str
#' @description Wrapper to fetch sequences from the archR result object as 
#' character
#'
#' @param res archR result object
#' 
#' @param iter Specify the iteration of archR result. If set to NULL 
#' (the default), the original set of sequences (`archRresult$rawSeqs`) is 
#' returned.
#' 
#' @param cl Specify the cluster number. Sequences belonging to this cluster in 
#' iteration `iter` of archR result are returned as character. When `iter` is 
#' NULL, this is treated as denoting the cluster number in archR's final
#' clustering solution (`archRresult$clustSol$clusters`). 
#' 
#' @param ord Specify TRUE if sequences are ordered by clusters. The original
#' ordering of the sequences can be fetched by setting `iter` to NULL and `ord` 
#' to FALSE.
#'
#' @details Setting iter to NULL will fetch sequences as per the final 
#' clustering solution of archR (`clustSol$clusters`). When `iter` is not 
#' NULL, use `cl` to further choose a particular cluster. When `cl` is NULL, 
#' the set of sequences returned can be ordered by clusters with `ord = TRUE`. 
#' Using `ord = FALSE` fetches the sequences by their original order. 
#' 
#' @return The selected DNA sequences from the DNAStringSet object as a 
#' character vector.
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' res <- system.file("extdata", "example_archRresult.rds", 
#'          package = "archR", mustWork = TRUE)
#'          
#' # Fetch sequences from 2nd cluster of archR's final solution
#' ans <- archR::seqs_str(res, iter=NULL, cl=2)
#' 
#' # Fetch all sequences ordered by the final clustering
#' ans <- archR::seqs_str(res, iter=NULL, cl=NULL, ord=TRUE)
#' 
#' # Fetch sequences belonging to first cluster in archR's first iteration
#' ans <- archR::seqs_str(res, iter=1, cl=1)
#' }
#' 
seqs_str <- function(res, iter = NULL, cl = NULL, ord = FALSE){
    
    ## when iter is NULL, return sequences belonging to a cluster (specified 
    ## by cl) in the final solution
    if(is.null(iter) && is.null(cl) && !ord){
        return(as.character(res$rawSeqs))
    }
    if(is.null(iter) && is.null(cl) && ord){
        use_ord <- unlist(res$clustSol$clusters)
        return(as.character(res$rawSeqs[use_ord]))
    }
    if(is.null(iter) && !is.null(cl)){
        cl_mem <- res$clustSol$clusters[[cl]]
        return(as.character(res$rawSeqs[cl_mem]))
    }
    if(!is.null(iter) && !is.null(cl)){
        clust_list <- get_seqs_clust_list(res$seqsClustLabels[[iter]])
        cl_mem <- clust_list[[cl]]
        return(as.character(res$rawSeqs[cl_mem]))
    }
    if(!is.null(iter) && is.null(cl) && ord){
        clust_list <- get_seqs_clust_list(res$seqsClustLabels[[iter]])
        use_ord <- unlist(clust_list)
        return(as.character(res$rawSeqs[use_ord]))
    }
    
}


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
## =============================================================================

.assign_samples_to_clusters <- function(clusterMembershipsVec, nClusters,
                                        iChunksColl, iChunkIdx) {
    returnClusterAsList <- vector("list", nClusters)
    for (i in seq_along(returnClusterAsList)) {
        returnClusterAsList[[i]] <-
            iChunksColl[[iChunkIdx]][clusterMembershipsVec == i]
    }
    return(returnClusterAsList)
}
## =============================================================================


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



## @title Collate clusters from hierarchical clustering
##
## @param clustList A list object giving indices of factors in each cluster
## @param globClustAssignments A list of cluster assignments
##
## @return collatedClustAssignments A list
.collate_clusters2 <- function(clustList, globClustAssignments,
                                flags = list(debugFlag = FALSE,
                                            verboseFlag = FALSE,
                                            plotVerboseFlag = FALSE,
                                            timeFlag = FALSE)) {
    if (is.null(clustList)) {
        ## When/if factor clustering is not performed, the 
        ## globClustAssingments variable is directly assigned
        collatedClustAssignments <- globClustAssignments
        return(collatedClustAssignments)
    }else{
        if(flags$debugFlag){
            message("=== This is what I got ===")
            message(.msg_print(globClustAssignments))
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
                message(.msg_print(clustList[[clustIdx]]))
                message("Size:", length(temp))
            }
            collatedClustAssignments[[clustIdx]] <- temp
        }
        if(flags$debugFlag){
            message("=== I am returning ===")
            message(.msg_print(collatedClustAssignments))
        }
        return(collatedClustAssignments)
    }
}
## =============================================================================

.msg_print <- function(vec_or_list){
    if(is.list(vec_or_list)){
        return(paste(vec_or_list, collapse = "\n"))
    }
    if(is.vector(vec_or_list)){
        return(paste(vec_or_list, collapse = ", "))
    }
}
## =============================================================================

.msg_pstr <- function(..., flg){
    if(flg){
        message(paste(...))
    }
}
## =============================================================================

## @title Update labels of sequences in a cluster
## @param oldSeqsClustLabels
##
## @param collatedClustAssignments
## @param flags List. Flags variable from archR config
##
## @return newSeqsClustLabels Vector os clustLabels
.update_cluster_labels <- function(oldSeqsClustLabels, 
                                collatedClustAssignments) {
    ## Important: The collatedClustAssignments variable is a list where each 
    ## element holds the indices of sequences falling in the same cluster, 
    ## like so: For a set of 200 sequences, collatedClustAssignments is
    ## [[1]]
    ## [1]   1   2   6  11  20  23  32  46  50  52  63  68  71  72  75  
    ## 80  82  83  85  86  87  88  92  99 102 104 105
    ##
    ## [[2]]
    ## [1]   3   5   7   8  10  12  13  16  17  18  19  22  25  30  31  
    ## 36  37  38  39  41  43  45  47  49  51  53  54
    ## --
    ## This is different than the globClustAssignment variable which will hold
    ## the indices of the factors that are combined into one cluster. 
    ## For instance, when there are 25 factors/clusters which have combined 
    ## to give 5 clusters, globClustAssignments will hold something like this:
    ## [[1]]
    ## [1] "3"  "13" "11" "16" "6" 
    ##
    ## [[2]]
    ## [1] "15" "1"  "8"  "25" "18"
    ## 
    .assert_archR_seqsClustLabels(oldSeqsClustLabels)
    .assert_archR_globClustAssignments(collatedClustAssignments)
    nClusters <- length(collatedClustAssignments)
    ## Use numerics w/ as.character and pre-sort to have them in the
    ## order that will be returned by levels (in get_seqs_clusters_in_list fn)
    candidateClustLabels <- sort(as.character(seq_len(nClusters)))
    newSeqsClustLabels <- oldSeqsClustLabels
    for (i in seq_len(nClusters)) {
        needUpdateIdx <- collatedClustAssignments[[i]]
        newSeqsClustLabels[needUpdateIdx] <-
            vapply(newSeqsClustLabels[needUpdateIdx], function(x) {
                        paste0(candidateClustLabels[i])}, character(1)
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
.prepare_chunks <- function(total_set, reqdChunkSize) {
    ## total_set is the set of seq_ids to be chunked (not array indices)
    if (is.null(total_set)) {
        stop("Preparing chunks, 'total_set' is NULL")
    }
    chunkLength <- length(total_set)
    if (chunkLength == 0) {
        stop("Preparing chunks, length of 'total_set' is 0")
    }
    .assert_archR_innerChunkSize_independent(reqdChunkSize)
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
        ##
        preparedChunks <- vector("list", length(chunkStarts))
        for (i in seq_along(chunkStarts)) {
            preparedChunks[[i]] <- 
                total_set[seq(chunkStarts[i],chunkEnds[i],by=1)]
        }
        ##
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
# 
# @param distMethod Specify the distance measure to be used. Default is cosine
# measure.
# 
# @param clustMethod Specify hc. This is the only option available.
# 
# @param linkage Specify one of linkage options for hclust.
# 
# @param flags Specify the flags from the config param set for archR.
# 
# @param returnOrder
# 
# @param useCutree By setting this arg to TRUE, the cutree version of 
# get_clusters func is used. This func uses silhouette and/or 
# Calinski-Harabasz index. See additional details there. Note: currently, this 
# is the only avaiable func so useCutree should beset to TRUE always, FALSE 
# will not work.
# 
# @param minClusters Passed to get_clusters function. See explanation there.
# 
# @param parentChunks
# 
# @return List with clustered factors
# 
# Change on 2020-12-28:
# - change default value for distMethod to 'euclid',
# - change default value for linkage to 'ward.D',
# - remove returnOrder arg. not needed; this was useful during development
# - 
# 
.handle_clustering_of_factors <- function(globFactorsMat,
                                        clustMethod = "hc",
                                        linkage = "average",
                                        distMethod = "cosangle",
                                        flags = list(debugFlag = FALSE,
                                                    verboseFlag = FALSE,
                                                    plotVerboseFlag = FALSE,
                                                    timeFlag = FALSE),
                                        useCutree = TRUE,
                                        minClusters = 2,
                                        parentChunks = NULL,
                                        ...) {
    ##
    .assert_archR_featuresMatrix(globFactorsMat)
    .assert_archR_flags(flags)
    dbg <- flags$debugFlag
    vrbs <- flags$verboseFlag
    ##
    gfDisMat <- .compute_factor_distances(globFactorsMat, distMethod)
    if(clustMethod == "hc"){
        as_dist_mat <- stats::as.dist(gfDisMat)
        temp_hclust <- stats::hclust(as_dist_mat, method = linkage)
        ord <- temp_hclust$order
        .msg_pstr("Cluster order by ", linkage," linkage w/ ",
            distMethod, " distance:", flg=dbg)
        .msg_pstr("New order: ", .msg_print(ord), flg=dbg)
        for(k in seq_len(length(ord)-1)){
            .msg_pstr(paste(paste(ord[k], ord[k+1], sep=", "), "=",
                    gfDisMat[ord[k], ord[k+1]], collapse="\n"), flg=dbg)
        }
        .msg_pstr("New order: ", .msg_print(ord), flg=dbg)
        .msg_pstr("=== Fetching clusters === ", flg=dbg)
        if(useCutree){
            clustList <- .get_clusters_from_hc_using_cutree(
                            hcObj = temp_hclust, distMat = as_dist_mat,
                            hStep = 0.05, parentChunks = parentChunks,
                            minClusters = minClusters, 
                            verbose = flags$debugFlag, ...)
        }
        
        if(length(clustList) == 1){
            clustList <- lapply(seq_len(max(clustList[[1]])), 
                                function(x){x})   
        }
        .msg_pstr("== DONE ==", "\nClustList: ", .msg_print(clustList), flg=dbg)
        return(clustList)
    }
}
## =============================================================================

.unfurl_nodeList <- function(nodeList){
    ##
    returnVal <- .assert_archR_list_properties(nodeList)
    if(returnVal != "FOO") stop(returnVal)
    ##
    element_lengths <- unlist(lapply(nodeList, function(elem){
        ifelse(is.list(elem), length(elem), 1)
        # else 1
    }))
    
    if(any(element_lengths != 1)){
        new_list <- vector("list", sum(element_lengths))
        iter1 <- 0
        iter2 <- 1
        while(iter2 <= length(nodeList)){
            if(!is.list(nodeList[[iter2]])){
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
## =============================================================================


############################## CUTREE VERSION ##################################
## This works best/is best referred with ward.D linkage
# @title Get clusters out of hierarchical clustering object using cutree. 
# Used internally by archR.
# 
# @description Clusters from a hierarchical clustering object are obtained 
# by using cutree at different heights of the tree. The optimum number of 
# clusters are decided based on the average silhouette value of the clustering.
# 
# @param hcObj The hierarchical clustering object as returned by 
# \code{\link[stats]{hclust}}.
# @param distMat The distance matrix that was used by hclust, a 
# \code{\link[stats]{dist}} object.
# @param hStep Numeric. The step size used to increment height values for 
# cutree. Default value is 0.05.
# @param parentChunks List. Specify the factor numbers in the previous 
# iteration of archR that factors in the current iteration resulted from. 
# Default value is NULL.
# @param keepSiblingsUncollated Logical. Specify TRUE if all clusters from a 
# single parent chunk should not be collated, even when recommended so based 
# on silhouette value computation. Default value FALSE. i.e. if they get 
# collated, it is left as collated.
# @param enableSwitchSilToCH Logical. When and if the clustering result has 
# singletons, using the average silhouette value may not be reliable to decide 
# on the optimum number of clusters. In this case, if this argument is TRUE, 
# the Calinski-Harabasz index is used instead. Setting it to FALSE, the
# switch is disabled, i.e., average silhouette value is used. Default is FALSE.
# @param minClusters Integer. Specify the minimum number of clusters to be 
# fetched from the HAC solution. 
# @param distThreshold Numeric. Specify a threshold of units of distance for 
# considering any two elements as close enough to be merged in the same 
# cluster. The default value is specified for Euclidean distance.
# @param verbose Logical. Specify TRUE for verbose output.
# 
# @importFrom stats cutree 
# @importFrom fpc calinhara
# @importFrom cluster silhouette
#
## Addtional comments:
## Keeping the min to 0 (as below) leads to as many clusters as 
## the number of elements. This errors when at mean[, "sil_width"]
## cut_heights <- seq(0, max(hcObj$height), by = hStep)
##
## Important: Address the question of whether any merging is required?
## -- We could use information from archR iteration
##    -- Obtain the cutree-using-h clustering here. 
##    -- Is it combining consecutive factors, like 1-2-3, 3-4 etc.? 
##          a. This could mean that it is undoing the clusters 
##          identified by archR in a previous iteration.
##          b. But, if 1 came from chunk 1 and 2-3 came from chunk 2, 
##          and if the clustering results in combining 1-2, then this 
##          is a plausible merge and should be left alone.
##          c. So, we use a parentChunks variable that notes the parent 
##          chunk for each of the current factors. This info can be 
##          used to make this decision unambiguously.
##          
## Update: 2020-12-13: We currently let these get combined if 
## HAC+cutree deems it fine, and we recommend that the clusters at the 
## iteration of archR can be left uncollated for reference. 
## The final stage will then hold the best possible set of clusters 
## guided by silhouette value. Now handled by logical argument 
## keepSiblingsUncollated.
##
## Update: 2020-12-24: With singletons in a clustering, using 
## silhouette values could be problematic, because it arbitrarily 
## assigns s(i) = 0 for a singleton cluster. This can lead to 
## smaller average silhouette value, artificially making that 
## clustering look bad. Thus, if and when we detect such a case, 
## test using the Calinski-Harabasz (CH) index automatically. 
## See argument `enableSwitchSilToCH`.
.get_clusters_from_hc_using_cutree <- function(hcObj, distMat, hStep = 0.05,
                                            parentChunks = NULL,
                                            keepSiblingsUncollated = FALSE,
                                            enableSwitchSilToCH = FALSE,
                                            minClusters = 2, 
                                            ## number of clusters in a 
                                            ## previous iteration on 
                                            ## which collation was performed
                                            distThreshold = 0.75,
                                            ## the default value for 
                                            ## - euclidean distance: 3
                                            ## - correlation distance: 0.75
                                            verbose = FALSE){
    ##
    if(minClusters < 2){
        warning("minClusters < 2. Setting it to 2")
        minClusters <- 2
    }
    .msg_pstr("minClusters is ", minClusters, flg=verbose)
    ##
    if(min(distMat) > distThreshold){
        clust_list <- lapply(hcObj$order, function(x) x)
        .msg_pstr("No element pairs close enough by given dist threshold: ", 
                distThreshold, flg=verbose)
        .msg_pstr("#Clusters: ", length(clust_list), flg=verbose)
        return(clust_list)
    }else{
        ##
        cut_heights <- seq(min(hcObj$height), max(hcObj$height), by = hStep)
        clust_list <- .get_clusts_sil_or_ch(cut_heights, hcObj, distMat, 
                        minClusters, use_sil = TRUE, verbose = verbose)
        ##
        .msg_pstr("#Clusts using sil.vals:", length(clust_list), flg=verbose)
        .msg_pstr(.msg_print(clust_list), flg=verbose)
        clust_list_lengths <- unlist(lapply(clust_list, length))
        if(enableSwitchSilToCH && any(clust_list_lengths == 1)){
            .msg_pstr("But singleton(s) present: clusters ", 
            paste(which(clust_list_lengths == 1), collapse= " "), flg=verbose)
            .msg_pstr("Using Calinski-Harabasz index instead", flg=verbose)
            ### Calinski-Harabasz Index
            clust_list <- .get_clusts_sil_or_ch(cut_heights, hcObj, distMat, 
                minClusters, use_sil = FALSE, verbose = verbose)
        }
        ## number of clusters decided, get final clustering result
        .msg_pstr("Final #Clusters: ", length(clust_list), flg=verbose)
        .msg_pstr(.msg_print(clust_list), flg=verbose)
        ##
        if(!is.null(parentChunks) && keepSiblingsUncollated){
            clust_list <- .check_and_uncollate_siblings(clust_list, 
                            parentChunks, verbose)
        }
        .msg_pstr("#Clusters: ", length(clust_list), flg=verbose)
        .msg_pstr(paste(clust_list), flg=verbose)
        return(clust_list)
    }
}
## =============================================================================

.check_and_uncollate_siblings <- function(clust_list, parentChunks, verbose){
    .msg_pstr("Uncollating siblings, in case...", flg=verbose)
    .msg_pstr("Checking parent chunks", flg=verbose)
    childrenPerParent <- lapply(unique(parentChunks), function(x){
        which(parentChunks == x)
    })
    updated_clust_list <- lapply(seq_along(clust_list), function(x){
        parents <- unique(parentChunks[clust_list[[x]]])
        thisX <- clust_list[[x]]
        if(length(parents) == 1 && length(thisX) > 1){
            ## if cluster elements come from same parent chunk of previous 
            ## iteration + now check if all elements in this cluster
            ## are the only children of that parent chunk (in other words 
            ## nothing got separated and combined with other factors)
            childrenThisParent <- childrenPerParent[[as.integer(parents)]]
            if(identical(as.numeric(childrenThisParent), as.numeric(thisX))){
                .msg_pstr("Parents are identical", flg=verbose)
                ## Update: 2020-12-13 
                ## See update above. 
                ## With keepSiblingsUncollated as FALSE, this check is not 
                ## performed.
                ## 
                ## Split the merge back into separate clusters
                return(lapply(childrenThisParent, function(y){y}))
            }else{
                .msg_pstr("returning thisX", flg=verbose)
                return(thisX)
            }
        }
        else{
            ## cluster elements come from different parent chunks 
            ## (of previous iteration)
            return(thisX)
        }
    })
    clust_list <- .unfurl_nodeList(updated_clust_list)
    clust_list
}
## =============================================================================

## This works best/is best referred with ward.D linkage
# @title Get clusters out of hierarchical clustering using either silhouette 
# value or Calinski-Harabasz index. 
# Used internally by archR.
# 
# @description Clusters from a hierarchical clustering object are obtained 
# by using cutree at different heights of the tree. The optimum number of 
# clusters are decided based on the average silhouette value of the clustering.
# 
# @param hcObj The hierarchical clustering object as returned by 
# \code{\link[stats]{hclust}}.
# @param distMat The distance matrix that was used by hclust, a 
# \code{\link[stats]{dist}} object.
# @param hStep Numeric. The step size used to increment height values for 
# cutree. Default value is 0.05.
# @param parentChunks List. Specify the factor numbers in the previous 
# iteration of archR that factors in the current iteration resulted from. 
# Default value is NULL.
# @param keepSiblingsUncollated Logical. Specify TRUE if all clusters from a 
# single parent chunk should not be collated, even when recommended so based 
# on silhouette value computation. Default value FALSE. i.e. if they get 
# collated, it is left as collated.
# @param enableSwitchSilToCH Logical. When and if the clustering result has 
# singletons, using the average silhouette value may not be reliable to decide 
# on the optimum number of clusters. In this case, if this argument is TRUE, 
# the Calinski-Harabasz index is used instead. Setting it to FALSE, the
# switch is disabled, i.e., average silhouette value is used. Default is FALSE.
# @param minClusters Integer. Specify the minimum number of clusters to be 
# fetched from the HAC solution. 
# @param distThreshold Numeric. Specify a threshold of units of distance for 
# considering any two elements as close enough to be merged in the same 
# cluster. The default value is specified for Euclidean distance.
# @param verbose Logical. Specify TRUE for verbose output.
# 
# @importFrom stats cutree 
# @importFrom fpc calinhara
# @importFrom cluster silhouette
#
.get_clusts_sil_or_ch <- function(cut_heights, hcObj, distMat, minClusters, 
                            use_sil = TRUE, verbose = FALSE){
    measure_cut_h <- unlist(lapply(cut_heights, function(x){
        foo_try <- stats::cutree(hcObj, h = x)
        names(foo_try) <- NULL
        if(length(unique(foo_try)) >= minClusters){
            if(use_sil){
                sils <- cluster::silhouette(foo_try, dist = distMat)
                retVal <- base::mean(sils[, "sil_width"])
            }else{
                retVal <- fpc::calinhara(distMat, foo_try)
            }
            retVal
        }else return(-100)
    }))
    
    ## multiple matches, first match index is returned with which.max
    cheight_idx <- which.max(measure_cut_h)
    if(use_sil){
        score_str <- "sil score"
    }else{
        score_str <- "CH index"
    }
    .msg_pstr("Max.", score_str, "at index =", cheight_idx, ", h:",
        cut_heights[cheight_idx], flg=verbose)
    ## number of clusters decided, get final clustering result
    cut_result <- stats::cutree(hcObj, h = cut_heights[cheight_idx])
    names(cut_result) <- NULL
    clust_list <- lapply(seq_along(unique(cut_result)),
        function(x) which(cut_result == x))
    clust_list
}
## =============================================================================


.regularizeMat <- function(basisMat, topN = 10){
    basisMat2 <- basisMat
    for(i in seq_len(ncol(basisMat))){
        asVec <- as.vector(basisMat[,i])
        threshold <- utils::tail(utils::head(
            sort(asVec, decreasing = TRUE), topN),1)
        basisMat2[(basisMat[,i] < threshold), i] <- 0.0
    }
    return(basisMat2)
}
## =============================================================================

#' @title Collate raw clusters at the chosen iteration of archR result
#' 
#' @description We use hierarchical clustering for reordering/collating raw
#' clusters from archR's given iteration.
#' 
#' @param result The archR result object.
#' 
#' @param iter Specify clusters at which iteration of archR are to be 
#' reordered/collated. Default is the last iteration of the archR result object.
#' 
#' @param clust_method Specify 'hc' for hierarchical clustering. Currently, only
#' hierarchical clustering is supported.
#' 
#' @param aggl_method One of linkage values as specified for hierarchical 
#' clustering with \code{\link[stats]{hclust}}. Default is 'ward.D'.
#' 
#' @param dist_method Distance measure to be used with hierarchical clustering. 
#' Available options are "euclid" (default), "cor" for correlation, "cosangle" 
#' for cosine angle, "modNW" for modified Needleman-Wunsch similarity (see 
#' \code{\link[TFBSTools]{PFMSimilarity}}). 
#' 
#' @param regularize Logical. Specify TRUE if regularization is to be performed 
#' before comparison. Default is FALSE. Also see argument 'topN'.
#' 
#' @param topn Use only the top N dimensions of each basis vector for 
#' comparing them. Note that since each basis vector has 4L or 16L (mono- or 
#' dinucleotides) dimensions, each dimension is a combination of nucleotide and 
#' its position in the sequence. This argument selects the top N dimensions of 
#' the basis vector. This is ignored when argument 'regularize' is FALSE.
#' 
#' @param collate Logical. Specify TRUE if collation using 
#' hierarchical agglomerative clustering is to be performed, otherwise FALSE. 
#' This argument is used by archR internally to obtain the ordering instead of 
#' collated clusters. 
#' 
#' @param flags Pass the flags object similar to the flags in configuration 
#' of the archR result object.
#'
#' @param ... ignored
#'
#' @return Returns the collated clusters. When `collate' is FALSE, 
#' it returns the already existing basis vectors, each as singleton clusters. 
#' The sequence cluster labels and sequence clusters are also handled 
#' accordingly.
#'  
#'
#' @importFrom stats hclust dist
#' @importFrom utils tail
#' @export
collate_archR_result <- function(result, 
                            iter = length(result$seqsClustLabels),
                            clust_method = "hc", aggl_method = "ward.D",
                            dist_method = "euclid", regularize = FALSE, 
                            topn = 50, collate = TRUE, flags, ...) {
    ##
    dbg <- flags$debugFlag
    vrbs <- flags$verboseFlag
    .msg_pstr("Collating clusters", flg=vrbs)
    stopifnot(iter > 0)
    basisMat <- get_clBasVec_m(result, iter=iter)
    clust_lab <- get_seqClLab(result, iter=iter)
    ##
    ## assert that topn value supplied is in the valid range
    if(!(topn > 0 && topn <= nrow(basisMat))){
        stop("Selected topn value ", topn," is outside expected range [", 
                .msg_print(c(1, nrow(basisMat))), "]")
    }
    if(!is.logical(regularize) && !is.logical(collate)){
        stop("Arguments 'regularize' and 'collate' expect a logical ",
                " value (TRUE/FALSE), found otherwise")
    }
    ##
    if(regularize) basisMat <- .regularizeMat(basisMat, topn)
    ##
    clust_lab_pre <- NULL
    if(iter > 1) clust_lab_pre <- get_seqClLab(result, iter=iter-1)
    parentChunks <- .get_parent_chunks(clust_lab, clust_lab_pre, iter, dbg)
    ##
    if(collate){
        factorsClustering <- 
            .handle_clustering_of_factors(basisMat, clustMethod = clust_method,
                    linkage = aggl_method, distMethod = dist_method,
                    flags = flags, parentChunks = parentChunks, ...)
        .msg_pstr("Factor clustering done, returning:", flg=dbg)
        .msg_pstr(paste("Returning: ", .msg_print(factorsClustering)), flg=dbg)
    }else{
        ## Do not reorder/collate, but prepare the return object
        nFactors <- get_clBasVec_k(result, iter) 
        factorsClustering <- vector("list", nFactors)
        factorsClustering <- lapply(seq_len(nFactors), function(x){x})
        .msg_pstr("No factor clustering, returning:", flg=dbg)
        .msg_pstr(.msg_print(factorsClustering), flg=dbg)
    }
    ##
    seqClusters <- get_seqs_clust_list(clust_lab)
    clusters <- .collate_clusters2(factorsClustering, seqClusters)
    clustLabels <- .update_cluster_labels(oldSeqsClustLabels = clust_lab, 
                                collatedClustAssignments = clusters)
    
    cluster_sol <- list(basisVectorsClust = factorsClustering,
                        clusters = clusters, seqsClustLabels = clustLabels)
    return(cluster_sol)
}
## =============================================================================

.get_parent_chunks <- function(clust_lab, clust_lab_pre, iter, flg){
    ## Prepare info on parent chunks for each cluster/factor at given iteration
    ## We get this info from the seqsClustLabels
    parentChunks <- NULL
    if(iter > 1){
        clustsThisIter <- sort(unique(clust_lab))
        parentChunks <- unlist(lapply(clustsThisIter, function(x){
            relSeqsIds <- which(clust_lab == x)
            unique(clust_lab_pre[relSeqsIds])
        }))
    }
    parentChunks
}
## =============================================================================


#' @title Retrieve sequence clusters as a list from the sequence labels
#' 
#' @description Given the sequence cluster labels from the archR result object,
#' returns the clusters separated as a list.
#'
#' @param seqs_clust_lab Sequences with cluster labels as in the archR result
#' object.
#'
#' @return A list holding sequence IDs belonging in each cluster.
#' 
#' @examples 
#' 
#' clustLabels <- sample(seq_len(4), 50, replace = TRUE)
#' print(clustLabels)
#' get_seqs_clust_list(clustLabels)
#' 
#' @export
get_seqs_clust_list <- function(seqs_clust_lab){
    ## check that labels are not empty/NULL
    .assert_archR_seqsClustLabels(seqs_clust_lab)
    clusterLevels <- levels(as.factor(seqs_clust_lab))

    seqs_clusters_as_a_list <- lapply(seq_along(clusterLevels),
                                        function(x){
                                            which(seqs_clust_lab ==
                                                    clusterLevels[x])
                                    })

    return(seqs_clusters_as_a_list)
}
## =============================================================================
