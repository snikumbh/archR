#' @title Handle directory creation
#'
#' @description Given the output directory name with its complete path, this
#' function checks if a directory of same name exists the given location. If
#' yes, then it adds a suffix (a number) to the given directory name, and
#' proceeds to create the directory.
#'
#' @param givenODir Specify the output directory name with its complete path.
#'
#'
#' @param flags List with four Logical elements as detailed.
#' \describe{
#'   \item{debugFlag}{Whether debug information for the run is printed}
#'   \item{verboseFlag}{Whether verbose information for the run is printed}
#'   \item{plotVerboseFlag}{Whether verbose plotting is performed for the run}
#'   \item{timeFlag}{Whether timing information is printed for the run}
#' }
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
## =============================================================================


## Getter function to fetch the features matrix from NMF result object
## (from python)
##Dependency on python script perform_nmf.py
get_features_matrix <- function(nmfResultObj){
    returnVal <- .assert_archR_list_properties(nmfResultObj)
    if (returnVal != "FOO") stop(returnVal)
    else return(as.matrix(nmfResultObj[[1]]))
}
## =============================================================================

## Getter function to fetch the samples matrix from NMF result object
## (from python)
## Dependency on python script perform_nmf.py
get_samples_matrix <- function(nmfResultObj){
    returnVal <- .assert_archR_list_properties(nmfResultObj)
    if (returnVal != "FOO") stop(returnVal)
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
# .get_hopach_cluster_medoidsIdx <- function(hopachObj){
#     .assert_archR_hopachObj(hopachObj, test_null = TRUE)
#     return(hopachObj$clustering$medoids)
# }
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
#' Default is 'stability'. Another option is 'cv' short for cross-validation.
#' Warning: The cross-validation approach can be time consuming and
#' computationally expensive than the stability-based approach. The CV approach
#' may be deprecate in the near future.
#' @param tol Numeric. Specify the tolerance value as criterion for choosing the
#' most appropriate number of NMF factors. Default is 1e-03.
#' @param bound Numeric. Specify the lower bound value as criterion for choosing
#' the most appropriate number of NMF factors. Default is 1e-08.
#' @param cvFolds Numeric. Specify the number of cross-validation folds used for
#'  model selection. Only used when modSelType is set to 'cv'.
#' @param parallelize Logical. Specify whether to parallelize the procedure.
#' Note that running archR serially can be time consuming. Consider
#' parallelizing with at least 2 or 4 cores. If Slurm is available, archR's
#' graphical user interface, accessed with \code{\link{run_archR_UI}}, enables 
#' provide all input data, set archR configuration, and directly
#' submit/monitor slurm jobs.
#' @param nCoresUse The number of cores to be used when `parallelize` is set
#' to TRUE. If `parallelize` is FALSE, nCoresUse is ignored.
#' @param nIterationsUse Numeric. Specify the number of bootstrapped iterations
#' to be performed with NMF.
#' @param alphaBase,alphaPow Specify the base value and the power for computing
#' 'alpha' in performing model selection for NMF. alpha = alphaBase^alphaPow.
#' Alpha specifies the regularization for NMF. Default: 0 and 1 respectively.
#' Warning: For future.
#' @param minSeqs Numeric. Specify the minimum number of sequences, such that
#' any cluster/chunk of size less than or equal to it will not be further
#' processed/clustered.
#' @param checkpointing Logical. Specify whether to write intermediate checkpoints
#' to disk as RDS files. Default is TRUE. Checkpoints and the final result are
#' saved to disk provided the oDir argument is set in \code{\link{archR}}.
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
                            kMax = 20,
                            modSelType = "stability",
                            tol = 10^-3,
                            bound = 10^-8,
                            cvFolds = 5,
                            parallelize = FALSE,
                            nCoresUse = NA,
                            nIterationsUse = 500,
                            alphaBase = 0,
                            alphaPow = 1,
                            minSeqs = 25,
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
    if(!requireNamespace("TFBSTools", quietly = TRUE)){
        warning("You chose 'modified Needleman-Wunsch' distance for 
            computing distances between NMF factors. This requires R 
            package 'TFBSTools', which is not available. Consider 
            installing the 'TFBSTools' package from Bioconductor and 
            re-run. Meanwhile, I am using the Euclidean distance.")
    }else{
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
}
## =============================================================================


## For hierarchical clustering object, return the cluster medoids
## We currently use the first element of the cluster as its medoid
.get_factors_from_factor_clustering2 <- function(listObj, globFactorsMat){
    ##
    .assert_archR_featuresMatrix(globFactorsMat)
    if (is.null(listObj)) {
        return(globFactorsMat)
    } else {
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
    #########################
    if(config$modSelType == "cv"){
        if(config$flags$debugFlag) {
            message("Performing cross validation-based model selection")
        }
        best_k <- .cv_model_select_pyNMF2(
                X = this_mat, param_ranges = config$paramRanges,
                kFolds = config$kFolds, parallelDo = config$parallelize,
                nCores = config$nCoresUse, nIterations = config$nIterationsUse,
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
        ## When order was changing:
        ## put samples matrix back in order it should be
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


intermediateResultsPlot <- function(seqsClustLabels, tss.seqs_raw = NULL,
                                    positions = NULL, iterVal = 0, fname = NULL){
## This function plots and prints resulting clusters -- the sequence image
## matrix (PNG file) and the sequence logos (PDF file).
## Input arguments: vector (size:nseqs) of cluster labels for sequence (for given iteration)
## Output: Nothing returned, files written to disk at specified location (fname)
##

if(is.numeric(iterVal)){
    name_suffix <- paste0("Iteration", iterVal)
    message("=== Intermediate Result ===")
}else{
    name_suffix <- paste0("Final")
    message("=== Final Result ===")
}
##
seqs_clusters_as_list_ordered <- get_seqs_clusters_in_a_list(seqsClustLabels)

message("Generating unannotated map of clustered sequences...")
image_fname <- paste0(fname, "ClusteringImage_", name_suffix, ".png")
message("Sequence clustering image written to: ", image_fname)
viz_seqs_as_acgt_mat_from_seqs(
    rawSeqs =  as.character(tss.seqs_raw[unlist(seqs_clusters_as_list_ordered)]),
                    position_labels = positions,
                    savefilename = image_fname,
                    fwidth = 450,
                    fheight = 900,
                    xt_freq = 5,
                    yt_freq = 100)

##
message("Generating architectures for clusters of sequences...")
arch_fname <- paste0(fname, "Architecture_SequenceLogos_", name_suffix, ".pdf")
message("Architectures written to: ", arch_fname)
plot_arch_for_clusters(
    tss.seqs_raw,
    list_of_elements = seqs_clusters_as_list_ordered,
    position_labels = positions,
    PDFfname = arch_fname)

}
## =============================================================================