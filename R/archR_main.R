#' @title
#' This function processes the given data set.
#'
#' @description Call this function to process a data set using archR.
#'
#' @param config archR configuration object as returned by
#' \code{\link{archRSetConfig}}. This is a required argument.
#' @param seqsMat A matrix of one-hot encoded sequences with sequences along
#' columns. This is a required argument.
#' @param seqsRaw A Biostrings::DNAStringSet object. The FASTA sequences as a
#' DNAStringSet object. This argument required argument.
#' @param seqsPositions Vector. Specify the tick labels for sequence positions.
#' Default is NULL.
#' @param thresholdItr Numeric. Specify the number of iterations to perform.
#' Default is 3.
#' @param setParsimony Logical vector. Specify if model selection by
#' cross-validation should prefer parsimonious solutions. Not required when
#' stability-based model selection is chosen. Length of the vector should match
#' number of threshold iterations specified in thresholdItr argument. TRUE
#' denotes parsimony is followed, FALSE otherwise.
#' @param setOCollation Logical vector. Specify for every iteration of archR
#' if collation of clusters from outer chunks should be performed. TRUE denotes
#' clusters are collated, FALSE otherwise.
#' @param fresh Logical. Specify if this is (not) a fresh run. Because
#' archR enables checkpointing, it is possible to perform additional iterations
#' upon clusters from an existing archR result (or a checkpoint) object. 
#' See UseOC argument. For example, when processing a set of FASTA sequences, 
#' if an earlier call to archR performed two iterations, and now you wish to 
#' perform a third, the arguments `fresh` and `UseOC` can be used. Simply set 
#' `fresh` to FALSE and assign the sequence clusters from iteration two from 
#' the earlier result to `UseOC`. As of v0.1.3, with this setting, archR 
#' returns a new result object as if the additional iteration performed is the 
#' only iteration.
#' @param UseOC List. Clusters to be further processed with archR. These can be 
#' from a previous archR result (in which case use archRresult$clustSol$clusters),
#' or simply clusters from any other method.
#' Warning: This has not been rigorously tested yet (v0.1.3).
#' @param oDir Character. Specify the output directory with its path. archR
#' will create this directory. If a directory with the given name exists at the
#' given location, archR will add a suffix to the directory name. This
#' change is reported to the user. Default is NULL. When NULL, just the result
#' is returned, and no plots or checkpoints or result is written to disk.
#'
#' @importFrom parallel makeCluster setDefaultCluster stopCluster
#'
#' @return A nested list of elements as follows:
#' \describe{
#' \item{seqsClustLabels}{A list with cluster labels for all sequences per 
#' iteration of archR. The cluster labels as stored as characters.}
#' 
#' \item{clustBasisVectors}{A list with information on NMF basis vectors per 
#' iteration of archR. Per iteration, there are two variables `nBasisVectors` 
#' storing the number of basis vectors after model selection,
#' and `basisVectors`, a matrix storing the basis vectors themselves. Dimensions
#'  of the `basisVectors` matrix are 4*L x nBasisVectors (mononucleotide 
#'  case) or 16*L x nBasisVectors (dinucleotide case).}
#'  
#' \item{clustSol}{The clustering solution obtained upon processing the raw 
#' clusters from the last iteration of archR's result. This is handled 
#' internally by the function \code{\link{reorder_archRresult}} using Euclidean 
#' distance and average linkage hierarchical clustering.}
#'  
#' \item{rawSeqs}{The input sequences as a DNAStringSet object.}
#' 
#' \item{timeInfo}{Stores the time taken (in minutes) for processing each 
#' iteration. This element is added only if `timeFlag` is set to TRUE in 
#' config.}
#' 
#' \item{config}{The configuration used for processing.}
#' \item{call}{The function call itself.}
#' } 
#' @export
archR <- function(config, seqsMat, seqsRaw, seqsPositions = NULL,
                  thresholdItr = 3,
                  setParsimony = c(FALSE, FALSE, FALSE),
                  setOCollation = c(TRUE, FALSE, FALSE),
                  fresh = TRUE,
                  UseOC = NULL,
                  oDir = NULL){
    ##
    ## assert thresholdItr is a positive integer
    if(!thresholdItr > 0) {
        stop("Expecting threshold iteration to be numeric and > 0")
    }
    ## better, provide a summary function
    ## message(utils::str(config))
    archRStartTime <- Sys.time()
    ##
    if(!is.null(oDir)){
        if(fresh){
            oDir <- handle_dir_creation(oDir, config$flags)
        }else{
            if(!dir.exists(oDir)){
                stop(oDir, " not found")
            }
        }
    }
    ##
    if(is.null(seqsPositions)){
        seqsPositions <- seq_len(Biostrings::width(seqsRaw[1]))
    }
    ##
    if(config$flags$plotVerboseFlag){
        allSequencesLogo <- plot_ggseqlogo_of_seqs(seqs = seqsRaw,
                               position_labels = seqsPositions,
                               title = paste("Sequence logo of all",
                                             length(seqsRaw),"sequences" ))
        ggsave(filename = file.path(oDir,"allSequencesLogo.pdf"),
               plot = allSequencesLogo,
               device = "pdf", width = 20, height = 2.5)
    }
    ## Make checks for params in configuration
    .assert_archR_config(config, ncol(seqsMat))
    .assert_archR_thresholdIteration(thresholdItr)
    ##
    ## ** To continue archR from an earlier run (further levels downstream)
    ## 1. Initializations of seqClustLabels, Factors etc should be
    ## appropriately handled.**
    ##
    ## Initialize sequence-cluster-labels and outerChunks
    seqsClustLabels <- rep("0", ncol(seqsMat))
    seqsClustLabelsList <- vector("list", thresholdItr)
    clustFactors <- vector("list", thresholdItr)
    architectures <- vector("list", thresholdItr)
    if(config$flags$timeFlag) timeInfo <- vector("list", thresholdItr)
    ##
    test_itr <- 1
    ##----------------------------------------------------------------------
    if(fresh){
        outerChunksColl <- vector("list", 1)
        ## Set outerChunks for first iteration
        outerChunksColl[[1]] <- seq(ncol(seqsMat))
        ##
    }else{
        message("Working on clusters from an earlier run")
        if(is.null(UseOC)){
            message("UseOC should not be NULL")
        }else{
            message("OK")
            ## UseOC is same as nxtOuterChunkColl
            seqsClustLabels <- .update_cluster_labels(seqsClustLabels,
                                                      UseOC)
            .assert_archR_OK_for_nextIteration(UseOC)
            outerChunksColl <- UseOC
            ##
        }
    }
    ##
    ##
    if(length(setParsimony) < thresholdItr){
        setParsimony <- rep(FALSE, thresholdItr)
        setParsimony[length(setParsimony)] <- TRUE
    }
    ##
    if(length(setOCollation) < thresholdItr){
        setOCollation <- rep(FALSE, thresholdItr)
        setOCollation[1] <- TRUE
        setOCollation[length(setParsimony)] <- FALSE #Earlier was TRUE
    }
    
    #### Start cluster only once
    if(config$parallelize){
        if(config$flags$debugFlag) message("Asked for parallelization")
        if(config$flags$debugFlag) message("Making cluster here")
        cl <- parallel::makeCluster(config$nCoresUse,
                                    type = "FORK")
        parallel::setDefaultCluster(cl)
        if(config$flags$verboseFlag) {
            message("Parallelization with ", config$nCoresUse , " cores")
        }
        if(config$flags$debugFlag) message(cl)
    }
    ####
    if(config$flags$verboseFlag) {
        if(config$modSelType == "cv"){
            message("Model selection by cross-validation")
        }
        if(config$modSelType == "stability"){
            message("Model selection by factor stability")
            message("Tolerance: ", config$tol, " & Bound: ", config$bound)
        }
    }

    ##
    ##--------------------------------------------------------------------------
    message("=== archR to perform ", thresholdItr, " iteration(s) ===")
    while (test_itr <= thresholdItr) {
        iterStartTime <- Sys.time()
        totOuterChunksColl <- length(outerChunksColl)
        if(config$flags$verboseFlag) {
            message("=== Iteration ", test_itr, ", ",
                totOuterChunksColl, " chunk(s) ===")
        }
        nxtOuterChunksColl <- vector("list")
        seqsClustLabels <- rep("0", ncol(seqsMat))
        intClustFactors <- NULL
        for (outerChunkIdx in seq_along(outerChunksColl)) {
            ##
            outerChunk <- outerChunksColl[[outerChunkIdx]]
            if(config$flags$verboseFlag) {
                message("[Outer chunk ", outerChunkIdx,"/",
                        totOuterChunksColl, "]",
                        " [Size: ", length(outerChunk), "]")
            }
            ## Make a decision to process based on size of chunk
            doNotProcess <- .decide_process_outer_chunk(config$minSeqs,
                                                        length(outerChunk),
                                                        config$kFolds)
            ##
            if (doNotProcess) {
                if(config$flags$verboseFlag) {
                    message("Decision: Skipping")
                }
                ## TO-DO: Could write function to manipulate
                ## collatedClustAssignments
                ## NOTE: Control enters here only from second iteration onwards,
                ## otherwise statements here can fail. Because, clustFactors is
                ## assumed to be populated already -- in the first iteration,
                ## test_itr is 0, and indexing using test_itr would fail.
                collatedClustAssignments <- list(outerChunk)
                ## We access clustFactors that are set in the previous iteration,
                ## when test_itr would be one less than it's current value.
                if (!is.null(intClustFactors)) {
                    intClustFactors <- cbind(intClustFactors,
                        as.matrix(
                        clustFactors[[test_itr-1]]$basisVectors[, outerChunkIdx])
                        )
                } else {
                    intClustFactors <- as.matrix(
                        clustFactors[[test_itr-1]]$basisVectors[, outerChunkIdx]
                        )
                }
            } else {
                if(config$flags$verboseFlag) {
                    message("Decision: Processing")
                }
                innerChunksColl <- .prepare_chunks(outerChunk,
                                                   config$innerChunkSize)
                ## Maintain these in a list, for collation later when
                ## all innerChunks in innerChunksColl have been processed
                globFactors <- vector("list", length(innerChunksColl))
                globClustAssignments <- vector("list", length(innerChunksColl))
                #################### INNER CHUNK FOR LOOP ######################
                for (innerChunkIdx in seq_along(innerChunksColl)) {
                    if (config$flags$verboseFlag) {
                        message("[Inner chunk ", innerChunkIdx, "/",
                            length(innerChunksColl), "]",
                            " [Size: ",
                            length(innerChunksColl[[innerChunkIdx]]),"]")
                    }
                    ##
                    ## Setting up sequences for the current chunk
                    this_seqsMat <-
                        seqsMat[, innerChunksColl[[innerChunkIdx]]]
                    if(config$modSelType == "cv" && !setParsimony[test_itr]) {
                        if(config$flags$verboseFlag)
                            message("=== Wihtout parsimony ===")
                    }
                    if(test_itr == 1 ||
                       length(outerChunk) > 0.9*config$innerChunkSize){
                        thisNMFResult <-
                            .handle_chunk_w_NMF2(innerChunkIdx,
                                        innerChunksColl,
                                        this_seqsMat,
                                        cgfglinear = TRUE,
                                        coarse_step = 10,
                                        monolinear = FALSE,
                                        askParsimony = setParsimony[test_itr],
                                        config)
                    }else{
                        ##
                        thisNMFResult <-
                            .handle_chunk_w_NMF2(innerChunkIdx,
                                        innerChunksColl,
                                        this_seqsMat,
                                        cgfglinear = TRUE,
                                        coarse_step = 5,
                                        monolinear = TRUE,
                                        askParsimony = setParsimony[test_itr],
                                        config)
                    }
                    .assert_archR_NMFresult(thisNMFResult)
                    globFactors[[innerChunkIdx]] <- thisNMFResult$forGlobFactors
                    globClustAssignments[[innerChunkIdx]] <-
                        thisNMFResult$forGlobClustAssignments
                } ## for loop over innerChunksColl ENDS here
                #################### INNER CHUNK FOR LOOP ######################
                ## We need globFactors, globClustAssignments
                ## Single unlist of globClustAssignments brings together
                ## clusters from different innerChunks into one collection
                globClustAssignments <- unlist(globClustAssignments,
                                                recursive = FALSE)
                ## CBind the factors from all inner chunks into one matrix
                globFactorsMat <- do.call(cbind, globFactors)
                ## These factors collected from all innerChunks may need
                ## clustering. This was earlier handled at the inner chunk
                ## level, but is now deferred to the outer chunk level.
                ## Therefore, the globFactorsClustering variable is directly
                ## set to NULL. Otherwise, this variable would hold a clustering
                ##  result object.
                globFactorsClustering <- NULL
                ## Manage factors
                if (!is.null(intClustFactors)) {
                    intClustFactors <-
                        cbind(intClustFactors,
                            .get_factors_from_factor_clustering2(
                                globFactorsClustering,
                                globFactorsMat))
                } else {
                    intClustFactors <-
                    .get_factors_from_factor_clustering2(globFactorsClustering,
                                                            globFactorsMat)
                }
                ## Manage collated cluster assignments
                collatedClustAssignments <-
                    .collate_clusters(globFactorsClustering,
                                        globClustAssignments)
            }  ## IfElse doNotProcess outer chunk ENDS
            ##
            .assert_archR_globClustAssignments(collatedClustAssignments)
            ## Assigning cluster labels can be done later, see below
            ## Collect (append) clusters at current level
            nxtOuterChunksColl <- append(nxtOuterChunksColl,
                                        collatedClustAssignments)
            ##
            chunksComplInfo <-
                paste0("[Outer chunk ", outerChunkIdx, "/", totOuterChunksColl,
                       "] complete")
            iterComplInfo <- paste0("[Iteration ", test_itr, "] complete")
            currInfo <- paste0("current total factors: ",
                               ncol(intClustFactors))
            nextIterInfo <- paste0("Current total chunks for next iteration: ",
                                length(nxtOuterChunksColl))
            ##
            if(config$flags$verboseFlag) {
                message(chunksComplInfo)
                if(outerChunkIdx == totOuterChunksColl)
                    message(iterComplInfo)
            }
            if(config$flags$debugFlag){
                # message(iterComplInfo)
                message(chunksComplInfo)
                message(currInfo)
                message(nextIterInfo)
            }
        }  ## for loop over outerChunksCollection ENDS
        if(config$flags$verboseFlag) {
            message("Managing clusters from outer chunk(s)")
        }

        ############### Managing clusters from outer chunks
        intClustFactorsClustering <- NULL

        ## Can intClustFactors ever be NULL?
        
        if(config$flags$debugFlag){
            message("Collation this iter:", setOCollation[test_itr],
                    "; #InnerChunks: ", length(innerChunksColl), 
                    "; #OuterChunks: ", totOuterChunksColl)
        }
        if (setOCollation[test_itr] && (length(innerChunksColl) > 1 
            || totOuterChunksColl > 1)) {
            ## ^The second condition in the IF statement protects against
            ## performing HAC/collation when #inner chunks & the number of 
            ## outer chunks is 1.
            if(config$flags$debugFlag) {
                message("Decision for outer chunk collation: Yes")
            }
            ## Cluster the factors using hierachical clustering
            ## Combinations considered:
            ## Cosine distance + Average Linkage
            ## Euclidean distance + Average Linkage
            ## modifiedNW distance + Average Linkage
            ## Euclidean distance + Complete Linkage
            ## Euclidean distance + Ward.D Linkage
            ## 
            ############## We use Complete linkage with Euclidean distance
            intClustFactorsClusteringEucCom <-
                .handle_clustering_of_factors(intClustFactors,
                                              clustMethod = "hc",
                                              linkage = "ward.D",## or "complete"?,
                                              distMethod = "euclid",
                                              flags = config$flags,
                                              returnOrder = FALSE)

            intClustFactors <- .get_factors_from_factor_clustering2(
                intClustFactorsClusteringEucCom, intClustFactors)
            ## Manage cluster assignments
            ## Call the hierarchical clustering version of collate_clusters
            ## This is named .collate_clusters2

            nxtOuterChunksColl <-
                .collate_clusters2(intClustFactorsClusteringEucCom,
                                   nxtOuterChunksColl, config$flags)

            ## Updating cluster labels for sequences can be done later after the
            ## clusters have been rearranged
            ## But we need this for rearrangement
            seqsClustLabels <- .update_cluster_labels(seqsClustLabels,
                                                      nxtOuterChunksColl,
                                                      flags = config$flags)
            ##
        } else {
            if(config$flags$debugFlag){
                message("Decision for outer chunk collation: No")
            }
            ## TO-DO: move this inside else block above?
            seqsClustLabels <- .update_cluster_labels(seqsClustLabels,
                                                      nxtOuterChunksColl,
                                                      flags = config$flags)
        }
        ############### MANAGING CLUSTERS FROM OUTER CHUNK ENDS ################
        seqsClustLabelsList[[test_itr]] <- seqsClustLabels
        ##
        clustFactors[[test_itr]] <-
            .setup_clustFactors_for_archR_result(intClustFactors)
        ##
        .assert_archR_OK_for_nextIteration(nxtOuterChunksColl)
        outerChunksColl <- nxtOuterChunksColl
        if(config$flags$plotVerboseFlag){
            if(!is.null(oDir)){
                intermediateResultsPlot(seqsClustLabelsList[[test_itr]],
                                    seqsRaw, positions = seqsPositions,
                                    iterVal = test_itr, fname = oDir)
            }
        }
        ##
        if(config$flags$verboseFlag) {
            itrComplMsg <- paste("Iteration", test_itr, "completed")
            itrComplTime <- ""
            if(config$flags$timeFlag){
                complTime <- format(as.numeric(Sys.time() - iterStartTime,
                                        units = "mins"), digits = 3)
                itrComplTime <- paste(" in", complTime, "mins")
                timeInfo[[test_itr]] <- as.double.difftime(complTime)
            }
            message(itrComplMsg, itrComplTime)
            archRComplMsg <- paste("ellapsed since start")
            archRComplTime <- ""
            if(config$flags$timeFlag){
                complTime <- format(as.numeric(Sys.time() - archRStartTime,
                                               units = "mins"), digits = 3)
                archRComplTime <- paste(complTime, "mins ")
                message(archRComplTime, archRComplMsg)
            }
        }
        ##
        ## test_itr was updated here
        ## write current iteration clustering to disk as RDS;
        ## good for checkpointing archR
        if(!is.null(oDir) && config$checkpointing && test_itr != thresholdItr){
            if(config$flags$verboseFlag) {
                message("Checkpointing set to TRUE. ",
                "Saving result of iteration ", test_itr, ".")
            }
            curr_archRresult <- list(seqsClustLabels = seqsClustLabelsList,
                                     clustBasisVectors = clustFactors,
                                     #architectures = architectures,
                                     rawSeqs = seqsRaw,
                                     config = config,
                                     call = match.call())
            rdsFilename <- paste0(oDir, "archRresult_checkpoint",
                                  test_itr, ".rds")
            saveRDS(curr_archRresult, file=rdsFilename)
        }
        ##
        test_itr <- test_itr + 1
    } ## algorithm while loop ENDS
    ##
    if(config$flags$timeFlag){
        temp_res <- list(seqsClustLabels = seqsClustLabelsList,
                        clustBasisVectors = clustFactors,
                        rawSeqs = seqsRaw,
                        timeInfo = timeInfo,
                        config = config,
                        call = match.call())
    }else{
        temp_res <- list(seqsClustLabels = seqsClustLabelsList,
                        clustBasisVectors = clustFactors,
                        rawSeqs = seqsRaw,
                        config = config,
                        call = match.call())
    }
    ##
    decisionToReorder <- TRUE
    iterations <- length(clustFactors)
    if(iterations > 1){
        prevIterNFactors <- clustFactors[[iterations-1]]$nBasisVectors
        currIterNFactors <- clustFactors[[iterations]]$nBasisVectors
        if(currIterNFactors == prevIterNFactors){
            decisionToReorder <- FALSE
            if(config$flags$debugFlag) message("Reordering decision: FALSE")
        }
    }
    ## For setting minClusters, note last iteration collated
    set_minClusters <- 2 ## the default value
    if(any(setOCollation)){
        lastItrC <- tail(which(setOCollation), 1)
        set_minClusters <- temp_res$clustBasisVectors[[lastItrC]]$nBasisVectors
    }
    ##
    temp_res_reord <- reorder_archRresult(temp_res,
                  iteration = thresholdItr,
                  clustMethod = "hc",
                  linkage = "ward.D", ## or average?
                  distMethod = "euclid",
                  minClusters = set_minClusters,
                  regularize = TRUE, 
                  topN = floor(0.5*length(seqsPositions)),
                  returnOrder = FALSE,
                  position_agnostic_dist = FALSE,
                  decisionToReorder = decisionToReorder,
                  config = temp_res$config)
    ## Print final stage output files to disk
    if(config$flags$plotVerboseFlag){
        if(!is.null(oDir)){
            intermediateResultsPlot(temp_res_reord$seqsClustLabels,
                                    seqsRaw, positions = seqsPositions,
                                    iterVal = "Final", fname = oDir)
        }
    }
    ##
    if(config$flags$timeFlag){
        temp_archRresult <- list(seqsClustLabels = seqsClustLabelsList,
                             clustBasisVectors = clustFactors,
                             clustSol = temp_res_reord,
                             rawSeqs = seqsRaw,
                             timeInfo = timeInfo,
                             config = config,
                             call = match.call())
    }else{
        temp_archRresult <- list(seqsClustLabels = seqsClustLabelsList,
                                 clustBasisVectors = clustFactors,
                                 clustSol = temp_res_reord,
                                 rawSeqs = seqsRaw,
                                 config = config,
                                 call = match.call())
    }
    .assert_archRresult(temp_archRresult)
    ##
    ## Write result to disk as RDS file
    if(!is.null(oDir)){
        rdsFilename <- paste0(oDir, "archRresult.rds")
        saveRDS(temp_archRresult, file=rdsFilename)
        ##
        message("=== Result saved to location ===")
        message(rdsFilename)
    }
    ##
    ## At the moment, archR result is not reordered, but this is left
    ## for later using reorder_archRresult function. This is because we
    ## currently use hierarchical clustering for this which can give different
    ## results with different distance measures. Performing this reordering
    ## later enables the user choose as per need.
    ##
    #### Stop cluster
    if(config$flags$debugFlag) {
        message("Stopping cluster...")
    }
    if(config$parallelize) parallel::stopCluster(cl)
    ####
    if(config$flags$verboseFlag) {
        archRComplMsg <- paste("archR exiting")
        archRComplTime <- ""
        if(config$flags$timeFlag){
            complTime <- format(as.numeric(Sys.time() - archRStartTime,
                                    units = "mins"), digits = 3)
            archRComplTime <- paste(",", complTime, "mins")
        }
        message(archRComplMsg, archRComplTime)
    }
    return(temp_archRresult)
}  ## archR function ENDS
## =============================================================================
