#' @title
#' This function processes the given data set.
#'
#' @description Call this function to process a data set using archR.
#'
#' @param config archR configuration object as returned by
#' \code{\link{archRSetConfig()}}. This is a required argument.
#' @param seqsMat A matrix of one-hot encoded sequences with sequences along
#' columns. This is a required argument.
#' @param seqRaw A Biostrings::DNAStringSet object. The FASTA sequences as a
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
#' using an existing archR result (or a checkpoint) object. See UseOC argument.
#' For example, when processing a set of FASTA sequences, if an earlier call
#' to archR performed two iterations, and now you wish to perform a third,
#' arguments fresh and UseOC can be used. Simply set fresh to FALSE and assign
#' the sequence clusters or iteration two from the earlier result to UseOC.
#' @param UseOC List. Clusters from an earlier iteration of archR or archR
#' result.
#' Warning: This has not been rigorously tested yet (v0.99.4).
#' @param oDir Character. Specify the output directory with its path. archR
#' will create this directory. If a directory with the given name exists at the
#' given location, archR will add a suffix to the directory name. This
#' change is reported to the user. Default is NULL. When NULL, just the result
#' is returned, and no plots or checkpoints or result is written to disk.
#'
#' @importFrom parallel makeCluster setDefaultCluster stopCluster
#'
#' @return Object of class archR (or at the moment, a list of lists)
#' @export
archR <- function(config, seqsMat, seqsRaw, seqsPositions = NULL,
                  thresholdItr = 3,
                  setParsimony = c(FALSE, FALSE, FALSE),
                  setOCollation = c(TRUE, FALSE, FALSE),
                  fresh = TRUE,
                  UseOC = NULL,
                  oDir = NULL){#"./archRresult") {
    ##
    ## assert thresholdItr is a positive integer
    stopifnot(thresholdItr > 0)
    ## better, provide a summary function
    str(config)
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
    ## ** To continue archR from an earlier run (further levels downstream)
    ## 1. Initializations of seqClustLabels, Factors etc should be
    ## appropriately handled.**
    ## Initialize sequence-cluster-labels and outerChunks
    seqsClustLabels <- rep("0", ncol(seqsMat))
    seqsClustLabelsList <- vector("list", thresholdItr)
    clustFactors <- vector("list", thresholdItr)
    architectures <- vector("list", thresholdItr)
    if(fresh){
        outerChunksColl <- vector("list", 1)
        ## Set outerChunks for first iteration
        outerChunksColl[[1]] <- seq(ncol(seqsMat))
    } else {
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
    ##
    test_itr <- 1
    ##--------------------------------------------------------------------------
    #### Start cluster only once
    if(config$parallelize){
        if(config$flags$debugFlag) message("Asked for parallelization")
        if(config$flags$debugFlag) message("Making cluster here")
        cl <- parallel::makeCluster(config$nCoresUse,
                                    type = "FORK",
                                    outfile = config$modSelLogFile)
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
    message("=== archR to perform ", thresholdItr, " iterations ===")
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
                ## clustering
                ############### HOPACH ON CLUSTERS FROM INNER CHUNK ############
                # msgSuffix1 <- ", Will skip collapsing of inner chunk clusters with HOPACH"
                # msgSuffix2 <- ", May collapse inner chunk clusters with HOPACH if decided"
                # if (test_itr == 1){
                # ## if (test_itr %% 2 != 1){
                #     message("Iteration: ", test_itr, msgSuffix1)
                #     hopachDecision <- FALSE
                # } else {
                #     if(config$flags$debugFlag) message("Iteration: ", test_itr, msgSuffix2)
                #     hopachDecision <- .decide_hopach(globFactorsMat,
                #                                 distMethod = "cosangle",
                #                                 withinMeasure = "mean")
                # }
                # hopachDecision <- FALSE
                #
                # if (hopachDecision) {
                #     ## Cluster the factors using hopach
                #     if(config$flags$debugFlag) message("Collapsing with HOPACH")
                #     globFactorsClustering <-
                #         .handle_clustering_of_factors(globFactorsMat,
                #                                     distMethod = "cosangle",
                #                                     flags = config$flags)
                # }
                ############### HOPACH ON CLUSTERS FROM INNER CHUNK ############
                globFactorsClustering <- NULL
                ## Manage factors
                if (!is.null(intClustFactors)) {
                    intClustFactors <-
                        cbind(intClustFactors,
                            .get_factors_from_factor_clustering(
                                globFactorsClustering,
                                globFactorsMat))
                } else {
                    intClustFactors <-
                    .get_factors_from_factor_clustering(globFactorsClustering,
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
        if (setOCollation[test_itr]) {
            if(config$flags$debugFlag) {
                message("SAMARTH: HOPACH yes for outer chunk collection")
            }
            ## Cluster the factors using hopach
            ## Combinations considered:
            ## Cosine distance + Average Linkage
            ## Euclidean distance + Average Linkage
            ## modifiedNW distance + Average Linkage
            ## Euclidean distance + Complete Linkage
            ##
            ############## We use Complete linkage with Euclidean distance
            intClustFactorsClusteringEucCom <-
                .handle_clustering_of_factors(intClustFactors,
                                              clustMethod = "hc",
                                              linkage = "complete",
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
            if(config$flags$debugFlag) message("Decision for outer chunk hopach: No")
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
    temp_archRresult <- list(seqsClustLabels = seqsClustLabelsList,
                        clustBasisVectors = clustFactors,
                        #architectures = architectures,
                        rawSeqs = seqsRaw,
                        config = config,
                        call = match.call())
    ## Write result to disk as RDS file
    if(!is.null(oDir)){
        rdsFilename <- paste0(oDir, "archRresult.rds")
        saveRDS(temp_archRresult, file=rdsFilename)
        ##
        message("=== Result saved to location ===")
        message(rdsFilename)
    }
    ##
    ## reordering returns the result object with an additional field
    ## message("Reordering archR clusters")
    #archRresult <- reorder_archRresult(temp_archRresult)
    ##
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
    #### Stop cluster
    if(config$flags$debugFlag) {
        message("Stopping cluster...")
    }
    if(config$parallelize) parallel::stopCluster(cl)
    ####
    return(temp_archRresult)
}  ## archR function ENDS
## =============================================================================
