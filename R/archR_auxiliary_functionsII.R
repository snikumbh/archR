#' @title Handle directory creation
#'
#' @export
handle_dir_creation <- function(givenODir, flags = list(debugFlag = FALSE,
                                                        verboseFlag = FALSE,
                                                        plotFlag = FALSE,
                                                        timeFlag = FALSE)){
    if(dir.exists(givenODir)){
        if(flags$verboseFlag) {
            message("-- Directory exists: -- ")
            message(givenODir)
            message("-- Changing name to: -- ")
        }
        allExistingDirs <- list.dirs(path = dirname(givenODir),
                                     recursive = FALSE)
        dirsThatMatch <- grep(pattern = basename(givenODir), allExistingDirs,
                              value = TRUE)
        ## Suffix an integer, because directory with given name (oDir)
        ## exists for length(dirsThatMatch) times
        name_suffix <- length(dirsThatMatch)
        givenODir <- paste(givenODir, name_suffix, sep = "_")
        while(dir.exists(givenODir)){
            name_suffix <- name_suffix + 1
            givenODir <- paste(givenODir, name_suffix, sep = "_")
        }
        if(flags$verboseFlag) {
            message(givenODir)
        }
    }
    retVal <- dir.create(paste0(givenODir, "/"), showWarnings = TRUE)
    stopifnot(retVal)
    returnODirName <- paste0(givenODir, "/")
    if(flags$verboseFlag) {
        message("-- Directory created for writing results -- ")
    }
    returnODirName
}



## Getter function to fetch the features matrix from NMF result object
## (from python)
##Dependency on python script perform_nmf.py
get_features_matrix <- function(nmfResultObj){
    returnVal <- .assert_archR_list_properties(nmfResultObj)
    if (returnVal != "SAMARTH") stop(returnVal)
    else return(as.matrix(nmfResultObj[[1]]))
}
## =============================================================================

## Getter function to fetch the samples matrix from NMF result object
## (from python)
## Dependency on python script perform_nmf.py
get_samples_matrix <- function(nmfResultObj){
    returnVal <- .assert_archR_list_properties(nmfResultObj)
    if (returnVal != "SAMARTH") stop(returnVal)
    else return(as.matrix(nmfResultObj[[2]]))
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
    .assert_archR_hopachObj(hopachObj, test_null = TRUE)
    return(hopachObj$clustering$medoids)
}
## =============================================================================

get_dimers_from_alphabet <- function(alphabet){
    if (!is.null(alphabet)) {
        return(do.call(paste0, expand.grid(alphabet, alphabet)))
    } else {
        stop("Expecting non-NULL alphabet")
    }
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
#' @param modSelType Character. Specify the model selection strategy to be used.
#' @param tol Numeric. Specify the tolerance value as criterion for choosing the
#' most appropriate number of NMF factors.
#' @param bound Numeric. Specify the lower bound value as criterion for choosing
#' the most appropriate number of NMF factors.
#' @param cvFolds Numeric. Specify the number of cross-validation folds used for
#'  model selection.
#' @param parallelize Logical. Specify whether to parallelize the procedure.
#' @param nCoresUse The number of cores to be used when `parallelize` is set to
#' TRUE. If `parallelize` is FALSE, nCoresUse is ignored.
#' @param nIterationsUse Specify the number of bootstrapped iterations to be
#' performed with NMF.
#' @param alphaBase,alphaPow Specify the base value and the power for computing
#' 'alpha' in performing model selection for NMF. alpha = alphaBase^alphaPow.
#' Alpha specifies the regularization for NMF. Default: 0 and 1 respectively.
#' @param minSeqs Specify the minimum number of sequences, such that any
#' cluster/chunk of size less than or equal to it will not further
#' processed/clustered.
#' @param modSelLogFile Specify a name for the file where model selection logs
#' will be wrtten.
#' @param checkpoint Logical. Specify whether to write intermediate checkpoints
#' to disk as RDS files. Default is TRUE.
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
                            modSelType = "stability", #other options "cv"
                            tol = 10^-4,
                            bound = 10^-8,
                            cvFolds = 5,
                            parallelize = FALSE,
                            nCoresUse = NA,
                            nIterationsUse = 500,
                            alphaBase = 0,
                            alphaPow = 1,
                            minSeqs = 25,
                            modSelLogFile = "log.txt",
                            checkpointing = TRUE,
                            flags = list(
                                debugFlag = FALSE,
                                timeFlag = FALSE,
                                verboseFlag = TRUE,
                                plotVerboseFlag = FALSE)
                            ) {
    ## Configuration Params that can be set by user
    archRconfig <- NULL
    archRconfig <- list(
                        modSelType = modSelType,
                        tol = tol,
                        bound = bound,
                        kFolds = cvFolds,
                        parallelize = parallelize,
                        nCoresUse = nCoresUse,
                        nIterationsUse = nIterationsUse,
                        paramRanges = list(
                            alphaBase = alphaBase,
                            alphaPow = alphaPow,
                            k_vals = seq(kMin, kMax, by = 1)
                        ),
                        innerChunkSize = innerChunkSize,
                        modSelLogFile = modSelLogFile,
                        checkpointing = checkpointing,
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
            stop("'minSeqs' should be at least 4 times 'kFolds'")
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
        message("Collating clusters")
        decision <- TRUE
    }
    if (!decision) message("No collation of clusters")
    return(decision)
}
## =============================================================================

## @param factorsMat A matrix holding the factors along the columns
## @param distMethod character A string specifying the distance measure to
## be computed. Values are: 'modNW' for modified Needleman Wunsch, and any
## distance measure that is possible with HOPACH
##
## Default value 'modNW'
##
## @return distance matrix from hopach (hdist object)
.compute_factor_distances <- function(factorsMat, distMethod = "modNW"){
    ## Assumption: Each column is a factor
    .assert_archR_featuresMatrix(factorsMat)
    if(distMethod == "modNW"){
        ## Turn the factors which are vectors into a 2D matrix of
        ## dinucs x positions
        dim_names <- get_dimers_from_alphabet(c("A", "C", "G", "T"))
        nPositions <- nrow(factorsMat)/length(dim_names)
        ##
        factorsMatList_as2D <- lapply(1:ncol(factorsMat),
                        function(x){matrix(factorsMat[,x],
                                       nrow = nrow(factorsMat)/nPositions,
                                       byrow = TRUE,
                                       dimnames = list(dim_names))
                        })
        ##
        factorsMatList_asPFMs <- lapply(1:length(factorsMatList_as2D),
            function(x){
                sinucSparse <- collapse_into_sinuc_matrix(
                                given_feature_mat = as.matrix(factorsMat[,x]),
                                dinuc_mat = factorsMatList_as2D[[x]],
                                feature_names = dim_names)
                sinucSparseInt <- matrix(as.integer(round(sinucSparse)),
                                nrow = 4, byrow = FALSE,
                                dimnames = list(rownames(sinucSparse)))
            })
        ##
        lenPFMs <- length(factorsMatList_asPFMs)
        scoresMat <- matrix(rep(0, lenPFMs*lenPFMs),
                                 nrow = lenPFMs)
        rownames(scoresMat) <- seq(1:nrow(scoresMat))
        colnames(scoresMat) <- seq(1:ncol(scoresMat))

        # relScoresMat <- scoresMat

        for(i in seq_len(lenPFMs)){
            for(j in seq_len(lenPFMs)){
                temp <- TFBSTools::PFMSimilarity(factorsMatList_asPFMs[[i]],
                                                 factorsMatList_asPFMs[[j]])
                scoresMat[i,j] <- temp["score"]
                # relScoresMat[i,j] <- temp["relScore"]
            }
        }
        ## currently we use scoresMat, so we only return that
        distMat <- max(scoresMat) - scoresMat
        return(distMat)
    } else{
        ## hopach::distancematrix function requires vectors along rows. Distances
        ## are computed between row vectors
        if (nrow(factorsMat) > ncol(factorsMat)) factorsMat <- t(factorsMat)
        hopachDistMat <- hopach::distancematrix(factorsMat, d = distMethod)
        ## hopachDistMat is a hopach hdist object
        stopifnot(hopachDistMat@Size == nrow(factorsMat))
        return(hopachDistMat)
    }

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


## for hierarchical clustering object
.get_factors_from_factor_clustering2 <- function(listObj, globFactorsMat){
    ##
    .assert_archR_featuresMatrix(globFactorsMat)
    if (is.null(listObj)) {
        return(globFactorsMat)
    } else {
        # .assert_archR_hopachObj(hopachObj, test_null = FALSE)
        # hopachMedoids <- .get_hopach_cluster_medoidsIdx(hopachObj)
        medoids <- unlist(lapply(listObj, function(x){x[1]}))
        return(as.matrix(globFactorsMat[ , medoids]))
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
.handle_chunk_w_NMF2 <- function(innerChunkIdx,
                                 innerChunksColl,
                                 this_mat,
                                 monolinear = FALSE,
                                 cgfglinear = TRUE,
                                 coarse_step = 10,
                                 askParsimony = TRUE,
                                 config){

    .assert_archR_flags(config$flags)
    ##

    if (is.null(this_mat) || !is.matrix(this_mat) &&
        !is(this_mat, "dgCMatrix")) {
        stop("Input matrix to model selection procedure is NULL or not a
            matrix")
    }
    ##
    if(config$modSelType == "cv"){
        if(config$flags$debugFlag) {
            message("Performing cross validation-based model selection")
        }
        best_k <- .cv_model_select_pyNMF2(
                X = this_mat, param_ranges = config$paramRanges,
                kFolds = config$kFolds, parallelDo = config$parallelize,
                nCores = config$nCoresUse, nIterations = config$nIterationsUse,
                logfile = config$modSelLogFile,
                set_verbose = 0, returnBestK = TRUE, #monolinear = monolinear,
                cgfglinear = cgfglinear, coarse_step = coarse_step,
                askParsimony = askParsimony
            )
    }
    #########################
    if(config$modSelType == "stability"){
        if(config$flags$debugFlag) {
            message("Performing stability-based model selection")
        }
        best_k <- .stability_model_select_pyNMF2(
            X = this_mat, param_ranges = config$paramRanges,
            parallelDo = config$parallelize, nCores = config$nCoresUse,
            nIterations = config$nIterationsUse,
            # logfile = config$modSelLogFile,
            tol = config$tol, bound = config$bound,
            flags = config$flags, returnBestK = TRUE, bootstrap = TRUE
        )
    }
    #########################
    if (best_k == max(config$paramRanges$k_vals)) {
        warning(c("WARNING: Best K for this subset == 'kMax'. ",
                  "Consider selecting a larger 'kMax' value, or\n",
                  "smaller innerChunkSize, or\n",
                  "perhaps, further increasing 'nIterationsUse'\n"),
                immediate. = TRUE)
    }
    if (config$flags$verboseFlag) {
        message("Best K for this subset: ", best_k)
    }

    ##
    if (config$flags$timeFlag) { start <- Sys.time() }
    ##
    if (best_k >= 1) {
        ## For fetching sequence clusters from samplesMat
        ## Cluster sequences
        ## New strategy, perform nRuns for bestK and use only the best one
        nRuns <- config$nIterationsUse
        if(config$flags$verboseFlag) {
            message("Fetching ", best_k, " clusters")
            # message("Fetching ", best_k, " clusters from best of ",
            #         nRuns," of NMF")
        }

        featuresMatrixList <- vector("list", nRuns)
        samplesMatrixList <- vector("list", nRuns)
        new_ord <- vector("list", nRuns)
        nmf_nRuns_list <- .perform_multiple_NMF_runs(X = this_mat,
                                             kVal = best_k,
                                             alphaVal = 0,
                                             parallelDo = config$parallelize,
                                             nCores = config$nCoresUse,
                                             nRuns = nRuns,
                                             bootstrap = TRUE)
        featuresMatrixList <- lapply(nmf_nRuns_list$nmf_result_list,
                                     function(x){
                                         get_features_matrix(x)
                                     })
        samplesMatrixList <- lapply(nmf_nRuns_list$nmf_result_list,
                                    function(x){
                                        get_samples_matrix(x)
                                    })
        ##
        new_ord <- nmf_nRuns_list$new_ord
        ## Get reconstruction accuracies for them
        bestQ2 <- -1
        for (nR in 1:nRuns){
            ##A <- this_mat[, new_ord[[nR]]]
            A <- this_mat
            recA <- as.matrix(featuresMatrixList[[nR]]) %*%
                as.matrix(samplesMatrixList[[nR]])
            this_q2 <- .compute_q2(as.matrix(A), recA)
            if(this_q2 > bestQ2){
                bestQ2 <- this_q2
                bestFeatMat <- featuresMatrixList[[nR]]
                bestSampMat <- samplesMatrixList[[nR]]
                bestOrd <- new_ord[[nR]]
            }
        }
        if (config$flags$debugFlag) {
            message("Best Q2 giving run found: ", bestQ2)
        }
        ##
        featuresMatrix <- bestFeatMat
        samplesMatrix <- bestSampMat
        # ## When order was chainging:
        # ## put samples matrix back in order it should be
        tempM <- bestSampMat
        samplesMatrix <- matrix(rep(NA, length(tempM)), nrow = nrow(tempM))
        samplesMatrix[ ,bestOrd] <- tempM
        #####
        if (config$flags$debugFlag) {
            message("Fetching ", best_k," cluster(s) w/ NMF scores")
        }
        clusterMembershipsForSamples <-
            .get_cluster_memberships_per_run(
                samplesMatrix = samplesMatrix,
                iChunksColl = innerChunksColl,
                iChunkIdx = innerChunkIdx)
        forGlobClustAssignments <- .assign_samples_to_clusters(
            clusterMembershipsVec =
                clusterMembershipsForSamples,
            nClusters = best_k,
            iChunkIdx = innerChunkIdx,
            iChunksColl = innerChunksColl)
        ###############
        #
    } else if (best_k < 1) {
        stop("Error in chosen number of factors: ", best_k)
    }
    .assert_archR_featuresMatrix(featuresMatrix)
    .assert_archR_globClustAssignments(forGlobClustAssignments)
    innerChunkNMFResult <- list(forGlobFactors = featuresMatrix,
                                forGlobClustAssignments =
                                    forGlobClustAssignments)
    ##
    return(innerChunkNMFResult)
}
## =============================================================================



.getMeanOfListOfMatrices <- function(listOfMats) {
    # Currently, assume all matrices have same dimensions
    # if discrepancy in dimensions, throws error "non-conformable arrays"
    meanMat <- Reduce("+", listOfMats)/length(listOfMats)
    return(meanMat)
}


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


reportArchFromResult <- function(archRresult, chooseLevel = NULL) {
    # INput archRresult
    # Output: architectures as sequence logos
    # Algorithm key points:
    # 1. Use cross_correlation and similarity measures from TFBSTools package
    # Stepwise algorithm:
    # 1. Get sequences in different clusters (at a given level/iteration; default: final iteration)
    # 2. Generate their sequence logos/PWMs
    # 3. Collect distinct PWMs as architectures in a archCollection object (list)
    # 4. Return the archCollection object (list of PWMs/sequence logos)
    if(is.null(chooseLevel)){
        # type is list, hence length
        chosenLevel <- length(archRresult$clustBasisVectors)
    }
    seqs_clusters_as_list_ordered <-
        get_seqs_clusters_in_a_list(archRresult$seqsClustLabels[[chosenLevel]])
    nClusters <- length(seqs_clusters_as_list_ordered)

    reportArch <- list()
    clusterPWMs <- vector("list", length(seqs_clusters_as_list_ordered))

    for(i in seq_along(seqs_clusters_as_list_ordered)){
        clusterPWMs[[i]] <- TFBSTools::toPWM(chen_tss.seqs_raw[seqs_clusters_as_list_ordered[[i]]],
                            type="prob")
    }
    # Fetch similarity scores of all pairs
    sims <- matrix(rep(0, nClusters*nClusters), nrow = nClusters, byrow = TRUE)
    for(i in seq_along(clusterPWMs)){
        for(j in i:length(clusterPWMs))
        sims[i,j] <- TFBSTools::PWMSimilarity(clusterPWMs[[i]], clusterPWMs[[j]],
                                   method="Euclidean")
    }



}


intermediateResultsPlot <- function(seqsClustLabels, tss.seqs_raw = NULL,
                                    positions = NULL, iterVal = 0, fname = NULL){
## This function plots and prints resulting clusters -- the sequence image
## matrix (PNG file) and the sequence logos (PDF file).
## Input arguments: vector (size:nseqs) of cluster labels for sequence (for given iteration)
## Output: Nothing returned, files written to disk at specified location (fname)
##
sorted_order <- sort(seqsClustLabels, index.return = TRUE)

##
seqs_clusters_as_list_ordered <- get_seqs_clusters_in_a_list(seqsClustLabels)
message("=== Intermediate Result ===")
message("Generating unannotated map of clustered sequences...")
image_fname <- paste0(fname, "ClusteringImage_Iteration", iterVal)
message("Sequence clustering image written to: ", image_fname, ".png")
viz_matrix_of_acgt_image(rawSeqs =  as.character(tss.seqs_raw[sorted_order$ix]),
                        position_labels = positions,
                        savefilename = image_fname,
                        fwidth = 450,
                        fheight = 900,
                        xt_freq = 5,
                        yt_freq = 100)

##
message("Generating architectures for clusters of sequences...")
arch_fname <- paste0(fname, "Architecture_SequenceLogos_Iteration", iterVal, ".pdf")
message("Architectures written to: ", arch_fname)
plot_arch_for_clusters_new(
    tss.seqs_raw,
    list_of_elements = seqs_clusters_as_list_ordered,
    position_labels = positions,
    PDFfname = arch_fname)

}
