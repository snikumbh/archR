get_features_matrix <- function(nmfResultObj){
    return(nmfResultObj[[1]])
}
## =============================================================================

get_samples_matrix <- function(nmfResultObj){
    return(nmfResultObj[[2]])
}
## =============================================================================

.get_hopach_cluster_medoidsIdx <- function(hopachObj){
    return(hopachObj$clustering$medoids)
}
## =============================================================================

get_dimers_from_alphabet <- function(alphabet){
    return(do.call(paste0, expand.grid(alphabet, alphabet)))
}
## =============================================================================

#' @title
#' Set archR run configuration in this function
#'
#' @description this will set the configuration
#'
#' @param innerChunkSize give
#' @param kMin give
#' @param kMax give
#' @param cvFolds give
#' @param parallelize give
#' @param nCoresUse give
#' @param nIterationsUse give
#' @param seedVal give
#' @param alphaBase give
#' @param alphaPow give
#' @param minSeqs give
#' @param modSelLogFile give
#' @param flags give
#'
#' @return a list with all params for archR set
#' @export
archRSetConfig <- function(innerChunkSize = 500,
                            kMin = 2,
                            kMax = 8,
                            cvFolds = 5,
                            parallelize = TRUE,
                            nCoresUse = 32,
                            nIterationsUse = 200,
                            seedVal = 10208090,
                            alphaBase = 0,
                            alphaPow = 1,
                            minSeqs = 25,
                            modSelLogFile = "log.txt",
                            flags = list(
                                debugFlag = FALSE,
                                timeFlag = FALSE,
                                verboseFlag = TRUE,
                                plotVerboseFlag = FALSE)
                            ) {
    ## Configuration Params that can be set by user
    archRconfig <- NULL
    archRconfig <- list(
                        kFolds = cvFolds,
                        parallelize = parallelize,
                        nCoresUse = nCoresUse,
                        nIterationsUse = nIterationsUse,
                        seedVal = seedVal,
                        paramRanges = list(
                            alphaBase = alphaBase,
                            alphaPow = alphaPow,
                            k_vals = seq(kMin, kMax, by = 1)
                        ),
                        innerChunkSize = innerChunkSize,
                        modSelLogFile = modSelLogFile,
                        minSeqs = minSeqs,
                        flags = list(
                            debugFlag = FALSE,
                            timeFlag = FALSE,
                            verboseFlag = TRUE,
                            plotVerboseFlag = FALSE
                        )
                    )
    return(archRconfig)
}
## =============================================================================

.decide_process_outer_chunk <-
    function(minThreshold, lengthOfOC, kFoldsVal) {
        # Assert that minThreshold > 4*kFoldsVal
        nFoldsCondition <- 4 * kFoldsVal
        base::stopifnot(minThreshold >= nFoldsCondition)
        doNotProcess <- FALSE
        if (lengthOfOC > 0) {
            if (lengthOfOC < minThreshold) {
                doNotProcess <- TRUE
                message("Sorry, will not process this small a chunk!")
            }
        } else {
            #doNotProcess <- TRUE
            message("WARNING: Outer chunk of size 0")
        }
        return(doNotProcess)
    }
## =============================================================================

.decide_hopach <- function(globFactorsMat,
                            distMethod = "cosangle",
                            withinMeasure = "mean") {
    ## Firstly:
    ## Very basically, if there are only two factors, we don't need HOPACH
    ## clustering of factors
    ## Secondly:
    ## If #factors > 2, we could need/do HOPACH, but if the distances between
    ## the factors are more or less similarly large (similar range), such that
    ## there is really no cluster/grouping, drop the idea of using HOPACH.
    globFactorsDistMat <- .compute_factor_distances(globFactorsMat,
                                                    distMethod = distMethod)
    ## Hopach suggestion: use same measure for 'within' and 'between'
    estClusters <- hopach::msscheck(globFactorsDistMat, within = "mean",
                            between = "mean")
    decision <- FALSE
    if (ncol(globFactorsMat) > 2 && estClusters[1] > 1) {
        message("Collating clusters from inner chunks")
        decision <- TRUE
    }
    if (!decision) message("No collation of clusters from inner chunks")
    return(decision)
}
## =============================================================================


.compute_factor_distances <- function(factorsMat, distMethod = "cosangle"){
    ## Ensure that entities to compare are along the rows, because
    ## Here, the factors are along the columns, and features along rows.
    ## hopach::distancematrix computes distances between rows of a matrix.
    ## Therefore, check and change if necessary.
    if (nrow(factorsMat) > ncol(factorsMat)) factorsMat <- t(factorsMat)
    distMat <- hopach::distancematrix(factorsMat, d = distMethod)
    ## distMat is a hopach hdist object
    assertthat::are_equal(distMat@Size, ncol(factorsMat))
    return(distMat)
}
## =============================================================================


.get_factors_from_factor_clustering <- function(hopachObj, globFactorsMat){
    ##
    if (is.null(hopachObj)) {
        return(globFactorsMat)
    } else {
        hopachMedoids <- .get_hopach_cluster_medoidsIdx(hopachObj)
        return(as.matrix(globFactorsMat[ , hopachMedoids]))
    }
}
## =============================================================================

.handle_chunk_w_NMF <- function(innerChunkIdx,
                                innerChunksColl,
                                this_tss.seqs,
                                globFactors,
                                globClustAssignments,
                                config) {
    ##
    ## On the given inner chunk,
    ## 1. perform model selection for #factors for NMF
    ## 2. Perform final NMF with chosen best_k (#Factors)
    ## 3. Store factors (globFactors)
    ## 4. Perform k-means clustering
    ## 5. Store cluster assignments (globClustAssignments)
    ## 6. Return updated globFactors, globClustAssignments
    ##
    if (config$flags$verboseFlag) {
        cat(paste0("Working on chunk: ", innerChunkIdx, " of ",
            length(innerChunksColl), " (chunkSize: ",
            ncol(this_tss.seqs), ") \n")
        )
    }
    ##
    model_selectK <-
        .cv_model_select_pyNMF(
            X = this_tss.seqs,
            param_ranges = config$paramRanges,
            kFolds = config$kFolds,
            parallelDo = config$parallelize,
            nCores = config$nCoresUse,
            nIterations = config$nIterationsUse,
            seed_val = config$seedVal,
            logfile = config$modSelLogFile,
            set_verbose = 0
        )
    if (config$flags$timeFlag) { print(Sys.time() - start) }
    ##
    best_k <- .get_best_K(model_selectK)
    ##
    if (config$flags$verboseFlag) {
        cat(paste0("Best K for this subset: ", best_k, "\n"))
    }
    if (config$flags$plotVerboseFlag) {
        q2_means_by_k_vals <-
            .get_q2_aggregates_chosen_var(model_selectK, model_selectK$k_vals,
                                            mean)
        Q2vsK <- .plot_cv_K(q2_means_by_k_vals)
        print(Q2vsK)
    }
    if (config$flags$verboseFlag ||
        config$flags$debugFlag) {
        cat("Performing NMF with K =", best_k, "\n")
    }
    ##
    if (config$flags$timeFlag) { start <- Sys.time() }
    ##
    result <- perform_nmf_func(
        this_tss.seqs,
        nPatterns = as.integer(best_k),
        nIter = as.integer(1000),
        givenAlpha = 0,
        givenL1_ratio = 1,
        seed_val = as.integer(config$seedVal)
    )
    if (config$flags$timeFlag) { print(Sys.time() - start) }
    ##
    featuresMatrix <- get_features_matrix(result)
    samplesMatrix <- get_samples_matrix(result)
    ##
    ## Collect factors for global list
    globFactors[[innerChunkIdx]] <- featuresMatrix
    ## Cluster sequences
    solKmeans <- get_clusters(samplesMatrix, clustMethod = "kmeans",
                                nCluster = best_k)
    if (config$flags$timeFlag) {print(Sys.time() - start)}
    ## Set the right cluster orders respective to the factor order in the chunk
    ## and, collect the clusters/cluster assignments for the global list
    globClustAssignments[[innerChunkIdx]] <-
        .map_clusters_to_factors(
            samplesMatrix = samplesMatrix,
            clustOrderIdx = solKmeans$reordering_idx,
            iChunksColl = innerChunksColl,
            iChunkIdx = innerChunkIdx,
            flags = config$flags
        )
    ##
    innerChunkNMFResult <- list(globFactors = globFactors,
                                globClustAssignments = globClustAssignments)
    ##
    return(innerChunkNMFResult)
}
## =============================================================================


assert_OK_for_nextIteration <- function(nxtOuterChunksColl) {
    ## Check any chunk of zero-length?
    if (any(lapply(nxtOuterChunksColl, length) == 0)) {
        message("WARNING: Chunks for next iteration have a problem")
        stop(which(lapply(nxtOuterChunksColl, length) == 0),
                " have zero lengths")
    }

}
