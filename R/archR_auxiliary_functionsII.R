## Getter function to fetch the features matrix from NMF result object
## (from python)
##Dependency on python script perform_nmf.py
get_features_matrix <- function(nmfResultObj){
    return(as.matrix(nmfResultObj[[1]]))
}
## =============================================================================

## Getter function to fetch the samples matrix from NMF result object
## (from python)
## Dependency on python script perform_nmf.py
get_samples_matrix <- function(nmfResultObj){
    return(as.matrix(nmfResultObj[[2]]))
}
## =============================================================================

## @title Get hopach cluster medoids
##
## @description Fetch the cluster medoids from hopach clustering result object
##
## @param hopachObj
##
## Dependency on hopach package
.get_hopach_cluster_medoidsIdx <- function(hopachObj){
    return(hopachObj$clustering$medoids)
}
## =============================================================================

get_dimers_from_alphabet <- function(alphabet){
    return(do.call(paste0, expand.grid(alphabet, alphabet)))
}
## =============================================================================

#' @title
#' Set archR run configuration
#'
#' @description This function sets the configuration for `archR`.
#'
#' @param innerChunkSize Numeric. Specify the size of the inner chunks of
#' sequences.
#' @param kMin Numeric. Specify the minimum of the range of values to be tested
#' for number of NMF basis vectors.
#' @param kMax Numeric. Specify the maximum of the range of values to be tested
#' for number of NMF basis vectors.
#' @param cvFolds Numeric. Specify the number of cross-validation folds used for
#'  model selection.
#' @param parallelize Logical. Specify whether to parallelize the procedure.
#' @param nCoresUse The number of cores to be used when `parallelize` is set to
#' TRUE. If `parallelize` is FALSE, nCoresUse is ignored.
#' @param nIterationsUse Specify the number of bootstrapped iterations to be
#' performed with NMF.
#' @param seedVal Specify the seed value for reproducibility.
#' @param alphaBase,alphaPow Specify the base value and the power for computing
#' 'alpha' in performing model selection for NMF. alpha = alphaBase^alphaPow.
#' Alpha specifies the regularization for NMF. Default: 0 and 1 respectively.
#' @param minSeqs Specify the minimum number of sequences, such that any
#' cluster/chunk of size less than or equal to it will not further
#' processed/clustered.
#' @param modSelLogFile Specify a name for the file where model selection logs
#' will be wrtten.
#' @param flags List with four Logical elements as detailed.
#' \describe{
#'   \item{debugFlag}{Whether debug information for the run is printed}
#'   \item{verboseFlag}{Whether verbose information for the run is printed}
#'   \item{plotVerboseFlag}{Whether verbose plotting is performed for the run}
#'   \item{timeFlag}{Whether timing information is printed for the run}
#' }
#'
#' @return a list with all params for archR set
#' @export
archRSetConfig <- function(innerChunkSize = 500,
                            kMin = 2,
                            kMax = 8,
                            cvFolds = 5,
                            parallelize = FALSE,
                            nCoresUse = NA,
                            nIterationsUse = 500,
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
                        flags = flags
                    )
    .assert_archR_config(archRconfig)
    return(archRconfig)
}
## =============================================================================


## @title Decide processing of outer chunk based on its size
##
## @description Function to make the decision on whether the given (outer)
## chunk should be processed
##
## @param minThreshold Numeric. This is the minSeqs param from archR config
## @param lengthOfOC Numeric. This is the length of the outer chunk
## @param kFoldsVal Numeric. This is the kFolds param in archR config
##
## If the lengthOfOC == 0, STOP
## If the minThreshold < 4*kFolds, STOP
## If the lengthOfOC < minThreshold, DO NOT PROCESS/return TRUE
##
.decide_process_outer_chunk <- function(minThreshold, lengthOfOC, kFoldsVal) {
        # Assert that minThreshold > 4*kFoldsVal
        nFoldsCondition <- 4 * kFoldsVal
        .assert_archR_minSeqs_independent(minThreshold)
        if (minThreshold < nFoldsCondition) {
            stop("'minSeqs' should be at least 4*'kFolds'")
        }
        # base::stopifnot(minThreshold >= nFoldsCondition)
        doNotProcess <- FALSE
        if (lengthOfOC > 0) {
            if (lengthOfOC < minThreshold) {
                doNotProcess <- TRUE
                message("Sorry, will not process this small a chunk: ",
                        lengthOfOC)
            }
        } else {
            # doNotProcess <- TRUE
            stop("Outer chunk of size 0")
        }
        return(doNotProcess)
    }
## =============================================================================

## This function makes the decision whether hopach-based processing should be
## performed (or will work)
##
## @param globFactorsMat A matrix object holding all factors as columns
## @param distMethod character A string specifying the method used for computing
## distance measure. Default 'cosangle'. Currently, 'cosangle' gives best
## results, so there are no other options here.
## @param withinMeasure character A string specifying whether "mean" or "median"
## should be used to compare between cluster dissimilarity. See more details in
## hopach.
##
## @return Logical TRUE/FALSE
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

## @param factorsMat A matrix holding the factors along the columns
## @param distMethod character A string specifying the distance measure to
## be computed.
## Default value 'cosangle'
##
## @return distance matrix from hopach (hdist object)
.compute_factor_distances <- function(factorsMat, distMethod = "cosangle"){
    ## Ensure that entities to compare are along the rows, because
    ## Here, the factors are along the columns, and features along rows.
    ## hopach::distancematrix computes distances between rows of a matrix.
    ## Therefore, check and change if necessary.
    .assert_archR_featuresMatrix(factorsMat)
    ## hopach::distancematrix function requires vectors along rows. Distances
    ## are computed between row vectors
    if (nrow(factorsMat) > ncol(factorsMat)) factorsMat <- t(factorsMat)
    distMat <- hopach::distancematrix(factorsMat, d = distMethod)
    ## distMat is a hopach hdist object
    stopifnot(distMat@Size == nrow(factorsMat))
    return(distMat)
}
## =============================================================================

# @title Get factors from factor clustering (hopach object)
#
# @description Returns the cluster medoid factors from hopach clustering of
# NMF factors
#
# @param hopachObj The hopach object holding hopach result
# @param globFactorsMat The global factors matrix
#
# @return If hopach object is not null, returns only the cluster medoid factors
# as a matrix, else, returns the complete factors matrix
.get_factors_from_factor_clustering <- function(hopachObj, globFactorsMat){
    ##
    .assert_archR_featuresMatrix(globFactorsMat)
    if (is.null(hopachObj)) {
        return(globFactorsMat)
    } else {
        .assert_archR_hopachObj(hopachObj, test_null = FALSE)
        hopachMedoids <- .get_hopach_cluster_medoidsIdx(hopachObj)
        return(as.matrix(globFactorsMat[ , hopachMedoids]))
    }
}
## =============================================================================

## @title Process a chunk wih NMF
##
##
## On the given inner chunk,
## 1. perform model selection for #factors for NMF
## 2. Perform final NMF with chosen best_k (#Factors)
## 3. Store factors (globFactors)
## 4. Fetch clusters using k-means clustering
##      - assign clusters <--> factors
## 5. Store cluster assignments (globClustAssignments)
## 6. Return updated globFactors, globClustAssignments
##
## - innerChunkIdx is needed to appropriately index into the global variables:
##    globFactors and globClustAssignments
## - globClustAssignments variable is updated inside the function to
## additionally hold new ones
## - Similarly, globFactors variable is updated inside the function to
## additionally hold new ones
.handle_chunk_w_NMF <- function(innerChunkIdx,
                                innerChunksColl,
                                this_mat,
                                config) {
    .assert_archR_flags(config$flags)
    ##
    if (config$flags$verboseFlag) {
        message("Working on inner chunk: ", innerChunkIdx, " of ",
                length(innerChunksColl), " [chunkSize: ", ncol(this_mat), "]")
    }
    ##
    if (is.null(this_mat) || !is.matrix(this_mat)) {
        stop("Input matrix to model selection procedure is NULL or not a
            matrix")
    }
    ############################### For fetching factors
    model_selectK <-
        .cv_model_select_pyNMF(
            X = this_mat, param_ranges = config$paramRanges,
            kFolds = config$kFolds, parallelDo = config$parallelize,
            nCores = config$nCoresUse, nIterations = config$nIterationsUse,
            seed_val = config$seedVal, logfile = config$modSelLogFile,
            set_verbose = 0
        )
    if (config$flags$timeFlag) { message(Sys.time() - start) }
    ##
    best_k <- .get_best_K(model_selectK)
    ##
    if (best_k == max(config$paramRanges$k_vals)) {
        warning(c("Best K for this subset == 'kMax'.",
                "Consider selecting a larger 'kMax' value, or\n",
                "smaller innerChunkSize, or\n",
                "perhaps, further increasing 'nIterationsUse'"),
                immediate. = TRUE)
    }
    if (config$flags$verboseFlag) {
        message("Best K for this subset: ", best_k)
    }
    if (config$flags$plotVerboseFlag) {
        q2_means_by_k_vals <- .get_q2_aggregates_chosen_var(
                            model_selectK, model_selectK$k_vals, mean)
        Q2vsK <- .plot_cv_K(q2_means_by_k_vals)
        print(Q2vsK)
    }
    # if (config$flags$verboseFlag || config$flags$debugFlag) {
    #     message("Performing NMF with K = ", best_k)
    # }
    ##
    if (config$flags$timeFlag) { start <- Sys.time() }
    ##
    result <- perform_nmf_func(this_mat, nPatterns = as.integer(best_k),
        nIter = as.integer(1000), givenAlpha = 0, givenL1_ratio = 1,
        seed_val = as.integer(config$seedVal)
        )
    if (config$flags$timeFlag) { print(Sys.time() - start) }
    ##
    featuresMatrix <- get_features_matrix(result)
    samplesMatrix <- get_samples_matrix(result)
    ##
    ## Will collect this as factors for global list
    ###############################
    ## For fetching sequence clusters from samplesMat
    ## Cluster sequences
    message("Fetching clusters")
    solKmeans <- get_clusters(samplesMatrix, clustMethod = "kmeans",
                                nCluster = best_k)
    if (config$flags$timeFlag) {print(Sys.time() - start)}
    ## Set the right cluster orders respective to the factor order in the chunk
    ## and, collect the clusters/cluster assignments for the global list
    forGlobClustAssignments <-
        .map_clusters_to_factors(
            samplesMatrix = samplesMatrix,
            clustOrderIdx = solKmeans$reordering_idx,
            iChunksColl = innerChunksColl, iChunkIdx = innerChunkIdx,
            flags = config$flags)
    ##
    .assert_archR_featuresMatrix(featuresMatrix)
    .assert_archR_globClustAssignments(forGlobClustAssignments)
    innerChunkNMFResult <- list(forGlobFactors = featuresMatrix,
                                forGlobClustAssignments =
                                    forGlobClustAssignments)
    ##
    return(innerChunkNMFResult)
}
## =============================================================================


## @title Setup the clustFactors list element for archR result object
##
## @description Function to set up the clustFactors variable for archR result
## object. Having a separate dedicated function enables seamless future changes.
##
## @param intClustFactors A matrix holding all factors from the just concluded
## iteration along the columns.
##
## @return A list with 2 elements having fixed element names. They are
## 'nBasisVectors':  This is the number of basis vectors (for easy info access)
##  and
## 'basisVectors': This is the matrix of basis vectors (intClustFactors itself)
.setup_clustFactors_for_archR_result <- function(intClustFactors) {
    returnList <- list(nBasisVectors = ncol(intClustFactors),
                        basisVectors = intClustFactors)
    return(returnList)
}
## =============================================================================
