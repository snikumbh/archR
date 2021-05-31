#' @title Handle directory creation
#'
#' @description Given the output directory name with its complete path, this
#' function checks if a directory of same name exists the given location. If
#' yes, then it adds a suffix (a number) to the given directory name, and
#' proceeds to create the directory.
#'
#' @param o_dir Specify the output directory name with its complete path.
#'
#' @param vrbs Set verbosity to TRUE or FALSE
#'
#' @return The (updated) dir name
#'
#' @export
handle_dir_creation <- function(o_dir, vrbs){
    ##
    cli::cli_alert_info(c("Output directory at path: ", 
                    "{.emph {dirname(o_dir)}}"))
    if(dir.exists(o_dir)){
        # .msg_pstr("-- Directory exists: -- ", o_dir, 
        #     "-- Changing name to: -- ", flg=vrbs)
        cli::cli_alert_warning("Directory exists: {.emph {basename(o_dir)}}")
        
        allExistingDirs <- list.dirs(path = dirname(o_dir),
                                        recursive = FALSE)
        dirsThatMatch <- grep(pattern = basename(o_dir), allExistingDirs,
                                value = TRUE)
        ## Suffix an integer, because directory with given name (oDir)
        ## exists for length(dirsThatMatch) times
        name_suffix <- length(dirsThatMatch)
        o_dir <- paste(o_dir, name_suffix, sep = "_")
        while(dir.exists(o_dir)){
            name_suffix <- name_suffix + 1
            o_dir <- paste(o_dir, name_suffix, sep = "_")
        }
        # .msg_pstr(o_dir, flg=vrbs)
        cli::cli_alert_info("Instead writing to: {.emph {basename(o_dir)}}")
    }
    retVal <- dir.create(paste0(o_dir, "/"), showWarnings = TRUE)
    stopifnot(retVal)
    returnODirName <- paste0(o_dir, "/")
    # .msg_pstr("-- Directory created for writing results -- ", flg=vrbs)
    # cli::cli_alert_success("")
    returnODirName
}
## =============================================================================


## Getter function to fetch the features matrix from NMF result object
## (from python)
##Dependency on python script perform_nmf.py
get_features_matrix <- function(nmfResultObj){
    returnVal <- .assert_archR_list_properties(nmfResultObj)
    if (returnVal != "FOO") stop(returnVal)
    return(as.matrix(nmfResultObj[[1]]))
}
## =============================================================================

## Getter function to fetch the samples matrix from NMF result object
## (from python)
## Dependency on python script perform_nmf.py
get_samples_matrix <- function(nmfResultObj){
    returnVal <- .assert_archR_list_properties(nmfResultObj)
    if (returnVal != "FOO") stop(returnVal)
    return(as.matrix(nmfResultObj[[2]]))
}
## =============================================================================

get_trimers_from_alphabet <- function(alph){
    if (is.null(alph)) stop("Expecting non-NULL alphabet")
    return(do.call(paste0, expand.grid(alph, alph, alph)))
    
}
## =============================================================================

get_dimers_from_alphabet <- function(alph){
    if (is.null(alph)) stop("Expecting non-NULL alphabet")
    return(do.call(paste0, expand.grid(alph, alph)))
    
}
## =============================================================================

manage_o_dir <- function(plt, o_dir){
    if(plt){
        if(is.null(o_dir)){
            stop(paste("'plot' flag is TRUE but 'o_dir' is not provided.",
                    "Did you forget to set 'o_dir'?"))
        }
    }
}
## =============================================================================

plot_all_seqs_logo <- function(seqs_raw, seqs_pos, dpath){
    if(is.null(dpath)){
        stop("directory path/name is NULL")
    }
    allSequencesLogo <- plot_ggseqlogo_of_seqs(
        seqs = seqs_raw,
        pos_lab = seqs_pos, 
        title = paste("Sequence logo of all", length(seqs_raw),"sequences" ))
    ##
    suppressWarnings(
        ggsave(filename = file.path(dpath, "allSequencesLogo.pdf"),
        plot = allSequencesLogo,
        device = "pdf", width = 20, height = 2.5))
}
## =============================================================================


#' @title
#' Set archR run configuration
#'
#' @description This function sets the configuration for `archR`.
#'
#' @param inner_chunk_size Numeric. Specify the size of the inner chunks of
#' sequences.
#' @param k_min Numeric. Specify the minimum of the range of values to be tested
#' for number of NMF basis vectors. Default is 1.
#' @param k_max Numeric. Specify the maximum of the range of values to be tested
#' for number of NMF basis vectors. Default is 20.
#' @param mod_sel_type Character. Specify the model selection strategy to 
#' be used. Default is 'stability'. Another option is 'cv', short for 
#' cross-validation. Warning: The cross-validation approach can be time 
#' consuming and computationally expensive than the stability-based approach. 
#' @param tol Numeric. Specify the tolerance value as criterion for choosing the
#' most appropriate number of NMF factors. Default is 1e-03. Current, this is
#' ignored.
#' @param bound Numeric. Specify the lower bound value as criterion for choosing
#' the most appropriate number of NMF factors. Default is 1e-08.
#' @param cv_folds Numeric. Specify the number of cross-validation folds used 
#' for model selection. Only used when mod_sel_type is set to 'cv'. Default 
#' value is 5.
#' @param parallelize Logical. Specify whether to parallelize the procedure.
#' Note that running archR serially can be time consuming. See `n_cores`. 
#' Consider parallelizing with at least 2 or 4 cores. If Slurm is available, 
#' archR's graphical user interface, accessed with \code{\link{run_archR_UI}},
#' enables providing all input data, setting archR configuration, and running 
#' archR directly by submitting/monitoring slurm jobs through the user 
#' interface.
#' @param n_cores The number of cores to be used when `parallelize` is set
#' to TRUE. If `parallelize` is FALSE, nCores is ignored.
#' @param n_iterations Numeric. Specify the number of bootstrapped iterations
#' to be performed with NMF. Default value is 100. When using cross-validation 
#' more than 100 (upto 500) iterations may be needed.  
#' @param alpha_base,alpha_pow Specify the base and the power for computing
#' 'alpha' in performing model selection for NMF. alpha = alpha_base^alpha_pow.
#' Alpha specifies the regularization for NMF. Default: 0 and 1 respectively.
#' _Warning_: Currently, not used (for future).
#' @param min_size Numeric. Specify the minimum number of sequences, such that
#' any cluster/chunk of size less than or equal to it will not be further
#' processed. Default is 25.
#' @param checkpointing Logical. Specify whether to write intermediate 
#' checkpoints to disk as RDS files. Checkpoints and the final result are 
#' saved to disk provided the `o_dir` argument is set in \code{\link{archR}}. 
#' When `o_dir` argument is not provided or NULL, this is ignored. 
#' Default is TRUE. 
#' @param flags List with four Logical elements as detailed.
#' \describe{
#'   \item{debug}{Whether debug information for the run is printed}
#'   \item{verbose}{Whether verbose information for the run is printed}
#'   \item{plot}{Whether verbose plotting is performed for the run}
#'   \item{time}{Whether timing information is printed for the run}
#' }
#'
#' @return a list with all params for archR set
#' 
#' @examples 
#' # Set archR configuration
#' archRconfig <- archR::archR_set_config(
#'     inner_chunk_size = 100,
#'     parallelize = TRUE,
#'     n_cores = 2,
#'     n_iterations = 100,
#'     k_min = 1,
#'     k_max = 20,
#'     mod_sel_type = "stability",
#'     tol = 10^-4,
#'     bound = 10^-8,
#'     flags = list(debug = FALSE, time = TRUE, verbose = TRUE,
#'         plot = FALSE)
#' )
#' 
#' 
#' @export
archR_set_config <- function(inner_chunk_size = 500,
                            k_min = 1,
                            k_max = 20,
                            mod_sel_type = "stability",
                            tol = 10^-3,
                            bound = 10^-8,
                            cv_folds = 5,
                            parallelize = FALSE,
                            n_cores = NA,
                            n_iterations = 500,
                            alpha_base = 0,
                            alpha_pow = 1,
                            min_size = 25,
                            checkpointing = TRUE,
                            flags = list(
                                debug = FALSE,
                                time = FALSE,
                                verbose = TRUE,
                                plot = FALSE)
                            ) {
    ## Configuration Params that can be set by user
    archRconfig <- NULL
    ##
    if(is.null(flags)){
        useFlags <- list(
            debugFlag = FALSE,
            timeFlag = FALSE,
            verboseFlag = TRUE,
            plotVerboseFlag = FALSE)
    }else{
        useFlags <- list(
            debugFlag = flags$debug,
            timeFlag = flags$time,
            verboseFlag = flags$verbose,
            plotVerboseFlag = flags$plot)
    }
    ##
    archRconfig <- list(modSelType = mod_sel_type,
                        tol = tol,
                        bound = bound,
                        kFolds = cv_folds,
                        parallelize = parallelize,
                        nCoresUse = n_cores,
                        nIterationsUse = n_iterations,
                        paramRanges = list(
                            alphaBase = alpha_base,
                            alphaPow = alpha_pow,
                            k_vals = seq(k_min, k_max, by = 1)
                        ),
                        innerChunkSize = inner_chunk_size,
                        checkpointing = checkpointing,
                        minSeqs = min_size,
                        flags = useFlags
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
        stopifnot(lengthOfOC > 0)
        stopifnot(minThreshold > 0)
        ## Assert that minThreshold > 4*kFoldsVal
        nFoldsCondition <- 4 * kFoldsVal
        .assert_archR_min_size_independent(minThreshold)
        if (minThreshold < nFoldsCondition) {
            stop("'min_size' should be at least 4 times 'kFolds'")
        }
        ##
        doNotProcess <- FALSE
        if (lengthOfOC <= minThreshold) {
            doNotProcess <- TRUE
            .msg_pstr("Sorry, will not process this small a chunk:",
                    lengthOfOC, flg=TRUE)
        }
        ##
        return(doNotProcess)
    }
## =============================================================================

## @param factorsMat A matrix holding the factors along the columns
## @param distMethod character A string specifying the distance measure to
## be computed. Values are: 'modNW' for modified Needleman Wunsch, and any
## distance measure that is possible with HOPACH
##
## Default value 'modNW'
## Edit on 2021-01-02:
## Default value changed to euclid instead of modNW which is computed using 
## a suggested package
##
## @return distance matrix from hopach (hdist object)
.compute_factor_distances <- function(factorsMat, distMethod = "euclid"){
    ## Assumption: Each column is a factor
    .assert_archR_featuresMatrix(factorsMat)
    ## Since the default distMethod is euclid/euclidean, when the user wishes
    ## to use any other distance methods/metrics, we can check if hopach exists
    ## If not, we ask the user to install it.
    ## This lets move hopach to Suggests
    distMethods_hopach <- c("cosangle", "abscosangle", 
                            "abseuclid", "cor", "abscor")
    distMethods_stats <- c("euclid")
    ##
    if(distMethod == "modNW"){
        distMat <- .get_modNW_dist(factorsMat)
        return(distMat)
    } 
    if(any(distMethod == distMethods_hopach)){
        distMat <- .get_hopach_dist(factorsMat, distMethod)
        return(distMat)
    }
    if(distMethod == distMethods_stats){
        distMat <- .get_stats_dist(factorsMat)
        return(distMat)
    }
}
## =============================================================================

.get_hopach_dist <- function(factorsMat, distMethod){
    if(!requireNamespace("hopach", quietly = TRUE)){
        stop("Please install R package 'hopach' to use ", distMethod, 
            " distance.")
    }else{
        ## - these distance metrics are available in hopach pkg
        ## - hopach::distancematrix func requires vectors along rows. 
        ## - Distances are computed between row vectors
        if (nrow(factorsMat) > ncol(factorsMat)){
            factorsMat <- t(factorsMat)
        }
        hopachDistMat <- hopach::distancematrix(factorsMat, d = distMethod)
        ## hopachDistMat is a hopach hdist object
        stopifnot(hopachDistMat@Size == nrow(factorsMat))
        ## make as.matrix as done for dist object in the 
        ## stats::dist case (see Else condition next)
        hopachDistMat <- hopach::as.matrix(hopachDistMat)
        return(hopachDistMat)
    }
}
## =============================================================================

.get_stats_dist <- function(factorsMat){
    ## dist method from stats // standard
    if (nrow(factorsMat) > ncol(factorsMat)) factorsMat <- t(factorsMat)
    as_dist <- stats::dist(factorsMat, method = "euclidean")
    distMat <- as.matrix(as_dist)
    return(distMat)
}

.get_modNW_dist <- function(factorsMat){
    if(!requireNamespace("TFBSTools", quietly = TRUE)){
        stop("Please install R package 'TFBSTools' for using modNW distance.")
    }else{
        ## Turn the factors which are vectors into a 2D matrix of
        ## dinucs x positions
        dim_names <- get_dimers_from_alphabet(c("A", "C", "G", "T"))
        nPositions <- nrow(factorsMat)/length(dim_names)
        ##
        # factorsMatList_as2D <- lapply(seq_len(ncol(factorsMat)),
        #     function(x){matrix(factorsMat[,x],
        #                     nrow = nrow(factorsMat)/nPositions,
        #                     byrow = TRUE,
        #                     dimnames = list(dim_names))
        #     })
        # ##
        # factorsMatList_asPFMs <- lapply(seq_len(length(factorsMatList_as2D)),
        #         function(x){
        #             sinucSparse <- collapse_into_sinuc_matrix(
        #                 given_feature_mat = as.matrix(factorsMat[,x]),
        #                 dinuc_mat = factorsMatList_as2D[[x]],
        #                 feature_names = dim_names)
        #             sinucSparseInt <- matrix(as.integer(round(sinucSparse)),
        #                 nrow = 4, byrow = FALSE,
        #                 dimnames = list(rownames(sinucSparse)))
        #         })
        ## After collapse_into_sinuc_matrix was updated, the above is changed
        ## to:
        factorsMatList_asPFMs <- lapply(seq_len(ncol(factorsMat)),
            function(x){
                ## A 16 x nPositions 2D matrix is obtained
                sinucSparse <- make_dinuc_PWMs(factorsMat[,x],
                                            add_pseudo_counts = FALSE,
                                            scale = FALSE)
                sinucSparseInt <- matrix(as.integer(round(sinucSparse)),
                                    nrow = 4, byrow = FALSE,
                                    dimnames = list(rownames(sinucSparse)))
            })
        
        ##
        lenPFMs <- length(factorsMatList_asPFMs)
        scoresMat <- matrix(rep(0, lenPFMs*lenPFMs), nrow = lenPFMs)
        rownames(scoresMat) <- seq(1,nrow(scoresMat),by=1)
        colnames(scoresMat) <- seq(1,ncol(scoresMat),by=1)
        
        for(i in seq_len(lenPFMs)){
            for(j in seq_len(lenPFMs)){
                temp <- TFBSTools::PFMSimilarity(
                    factorsMatList_asPFMs[[i]], factorsMatList_asPFMs[[j]])
                scoresMat[i,j] <- temp["score"]
                # relScoresMat[i,j] <- temp["relScore"]
            }
        }
        ## currently we use scoresMat, so we only return that
        distMat <- max(scoresMat) - scoresMat
    }
    return(distMat)
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
## 
## doRegularize arg is passed on to stability_model _selection function
## For first iteration, doRegularize will be FALSE, and TRUE in subsequent 
## iterations of archR
## 
.handle_chunk_w_NMF2 <- function(innerChunkIdx,
                                    innerChunksColl,
                                    this_mat,
                                    monolinear = FALSE,
                                    cgfglinear = TRUE,
                                    coarse_step = 10,
                                    askParsimony = TRUE,
                                    doRegularize = FALSE,
                                    config, oDir, test_itr, oChunkIdx){
    .assert_archR_flags(config$flags)
    dbg <- config$flags$debugFlag
    vrbs <- config$flags$verboseFlag
    tym <- config$flags$timeFlag
    ##
    if (is.null(this_mat) || !is.matrix(this_mat) &&
        !is(this_mat, "dgCMatrix")) {
        stop("Input matrix to model selection procedure is NULL or not a
            matrix")
    }
    #########################
    if(config$modSelType == "cv"){
        suffix_str <- ", without parsimony"
        if(askParsimony) suffix_str <- ", with parsimony"
        .msg_pstr("Performing cross validation-based model selection",
                suffix_str, flg=dbg)
        best_k <- .cv_model_select_pyNMF2(
                X = this_mat, param_ranges = config$paramRanges,
                kFolds = config$kFolds, parallelDo = config$parallelize,
                nCores = config$nCoresUse, nIterations = config$nIterationsUse,
                verboseFlag = config$flags$verboseFlag,
                debugFlag = config$flags$debugFlag,
                returnBestK = TRUE, cgfglinear = cgfglinear, 
                coarse_step = coarse_step,
                askParsimony = askParsimony
            )
    }
    #########################
    if(config$modSelType == "stability"){
        .msg_pstr("Performing stability-based model selection", flg=dbg)
        if(doRegularize) .msg_pstr("Regularization w/ dispersion", flg=dbg)
        best_k <- .stability_model_select_pyNMF2(
            X = this_mat, param_ranges = config$paramRanges,
            parallelDo = config$parallelize, nCores = config$nCoresUse,
            nIterations = config$nIterationsUse, tol = config$tol, 
            bound = config$bound, flags = config$flags, 
            returnBestK = TRUE, bootstrap = TRUE
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
    .msg_pstr("Best K for this subset:", best_k, flg=vrbs)
    
    
    ##
    if (best_k >= 1) {
        ## For fetching sequence clusters from samplesMat
        ## Cluster sequences
        ## New strategy, perform nRuns for bestK and use only the best one
        nRuns <- config$nIterationsUse
        # .msg_pstr("Fetching ", best_k, " clusters", flg=(vrbs || dbg))

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
                                        get_features_matrix)
        samplesMatrixList <- lapply(nmf_nRuns_list$nmf_result_list,
                                        get_samples_matrix)
        ##
        new_ord <- nmf_nRuns_list$new_ord
        ## Get reconstruction accuracies for them
        bestQ2 <- -1
        for (nR in seq_len(nRuns)){
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
        # .msg_pstr("Best Q2 giving run found: ", bestQ2, flg=dbg)
        # cli::cli_alert_info("Fetched {best_k} clusters")
        ##
        featuresMatrix <- bestFeatMat
        samplesMatrix <- bestSampMat
        ## When order was changing:
        ## put samples matrix back in order it should be
        tempM <- bestSampMat
        samplesMatrix <- matrix(rep(NA, length(tempM)), nrow = nrow(tempM))
        samplesMatrix[ ,bestOrd] <- tempM
        #####
        # .msg_pstr("Fetching ", best_k," cluster(s)", flg=dbg)
        clusterMembershipsForSamples <-
            .get_cluster_memberships_per_run(samplesMatrix = samplesMatrix,
                iChunksColl = innerChunksColl, iChunkIdx = innerChunkIdx, 
                oDir, test_itr, oChunkIdx)
        ##
        ## Could handle overfitting here
        if(best_k > 1){
            has_overfit <- .detect_overfitting(samplesMatrix, 
                                                clusterMembershipsForSamples, 
                                                minSeqs = 50)
            # print("Overfit at:")
            # if(length(has_overfit) > 0){ 
            #     print(paste("Overfit at:", has_overfit))
            # }else{
            #     print("No Overfit")
            # }
            ## 
            ## -- Note which clusters are overfit
            ## -- Remove those columns from featuresMat
            ## -- Adjust clustMemberships
            if(length(has_overfit) > 0){
                # print(has_overfit)
                clusterMembershipsForSamples <- 
                    .adjustSampleMemberships(clusterMembershipsForSamples, 
                                        samplesMatrix, has_overfit)
                featuresMatrix <- as.matrix(featuresMatrix[, -c(has_overfit)])
                best_k <- best_k - length(has_overfit)
                # .msg_pstr("Adjusting for overfitting, fetched ", 
                #     best_k, "cluster(s)", flg=vrbs)
                cli::cli_alert_info(c("Adjusting for overfitting, ",
                                    "fetched {best_k} cluster{?s}"))
            }
            # else{
            #     print("None")
            # }
        }
        ##
        forGlobClustAssignments <- .assign_samples_to_clusters(
            clusterMembershipsVec = clusterMembershipsForSamples,
            nClusters = best_k, iChunkIdx = innerChunkIdx, 
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


intermediateResultsPlot <- function(seq_lab, seqs_raw = NULL,
                                pos_lab = NULL, iter = 0, fname = NULL,
                                name_suffix = NULL,
                                vrbs = TRUE){
## This function plots and prints resulting clusters -- the sequence image
## matrix (PNG file) and the sequence logos (PDF file).

if(is.null(name_suffix)){
    if(is.numeric(iter)){
        name_suffix <- paste0("Iteration", iter)
        # .msg_pstr("=== Intermediate Result ===", flg=vrbs)
        cli::cli_rule(left="Intermediate result")
    }else{
        name_suffix <- paste0("Final")
        # .msg_pstr("=== Final Result ===", flg=vrbs)
        cli::cli_rule(left="Final result")
    }
}else{
    .msg_pstr("=== On-demand Result ===", flg=vrbs)
}
##
cli::cli_alert_info("Output directory: {.emph {fname}}")
    
seqs_clust_list_ord <- get_seqs_clust_list(seq_lab)
seqs_clust_vec_ord <- unlist(seqs_clust_list_ord)
# .msg_pstr("Generating unannotated map of clustered sequences...", flg=vrbs)
image_fname <- paste0(fname, "ClusteringImage_", name_suffix, ".png")
# .msg_pstr("Sequence clustering image written to: ", image_fname, flg=vrbs)
cli::cli_alert_info(c("Sequence clustering image written to: ", 
                        "{.emph {basename(image_fname)}}"))
viz_seqs_acgt_mat_from_seqs(
    seqs =  as.character(seqs_raw[seqs_clust_vec_ord]),
                    pos_lab = pos_lab,
                    save_fname = image_fname,
                    f_width = 450,
                    f_height = 900,
                    xt_freq = 5,
                    yt_freq = 100)

##
# .msg_pstr("Generating architectures for clusters of sequences...", flg=vrbs)

arch_fname <- paste0(fname, "Architecture_SequenceLogos_", name_suffix, ".pdf")
# .msg_pstr("Architectures written to: ", arch_fname, flg=vrbs)
cli::cli_alert_info(c("Architectures written to: ", 
                        "{.emph {basename(arch_fname)}}"))

plot_arch_for_clusters(seqs = seqs_raw, 
                        clust_list = seqs_clust_list_ord,
                        pos_lab = pos_lab,
                        pdf_name = arch_fname)
}
## =============================================================================