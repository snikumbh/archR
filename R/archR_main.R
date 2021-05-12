#' @title
#' This function processes the given data set.
#'
#' @description Call this function to process a data set using archR.
#'
#' @param config archR configuration object as returned by
#' \code{\link{archR_set_config}}. This is a required argument.
#' @param seqs_ohe_mat A matrix of one-hot encoded sequences with sequences 
#' along columns. This is a required argument.
#' @param seqs_raw A \code{\link[Biostrings]{DNAStringSet}} object. The FASTA 
#' sequences as a DNAStringSet object. This argument required argument.
#' @param seqs_pos Vector. Specify the tick labels for sequence positions.
#' Default is NULL.
#' @param threshold_itr Numeric. Specify the number of iterations to perform.
#' Default is 3.
#' @param set_parsimony Logical vector. Specify if model selection by
#' cross-validation should prefer parsimonious solutions. Not required when
#' stability-based model selection is chosen. Length of the vector should match
#' number of iterations specified in 'threshold_itr' argument. TRUE denotes 
#' parsimony is followed, FALSE otherwise.
#' @param set_ocollation Logical vector. Specify for every iteration of archR
#' if collation of clusters from outer chunks should be performed. TRUE denotes
#' clusters are collated, FALSE otherwise.
#' @param fresh Logical. Specify if this is (not) a fresh run. Because
#' archR enables checkpointing, it is possible to perform additional iterations
#' upon clusters from an existing archR result (or a checkpoint) object. 
#' See 'use_oc' argument. For example, when processing a set of FASTA sequences,
#' if an earlier call to archR performed two iterations, and now you wish to 
#' perform a third, the arguments `fresh` and `use_oc` can be used. Simply set 
#' `fresh` to FALSE and assign the sequence clusters from iteration two from 
#' the earlier result to `use_oc`. As of v0.1.3, with this setting, archR 
#' returns a new result object as if the additional iteration performed is the 
#' only iteration.
#' @param use_oc List. Clusters to be further processed with archR. These can 
#' be from a previous archR result (in which case use 
#' archRresult$clustSol$clusters), or simply clusters from any other method.
#' Warning: This has not been rigorously tested yet (v0.1.3).
#' @param o_dir Character. Specify the output directory with its path. archR
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
#' and `basisVectors`, a matrix storing the basis vectors themselves. 
#' Dimensions of the `basisVectors` matrix are 4*L x nBasisVectors 
#' (mononucleotide case) or 16*L x nBasisVectors (dinucleotide case).}
#'  
#' \item{clustSol}{The clustering solution obtained upon processing the raw 
#' clusters from the last iteration of archR's result. This is handled 
#' internally by the function \code{\link{collate_archR_result}} using the
#' default setting of Euclidean distance and ward.D linkage hierarchical 
#' clustering.}
#'  
#' \item{rawSeqs}{The input sequences as a DNAStringSet object.}
#' 
#' \item{timeInfo}{Stores the time taken (in minutes) for processing each 
#' iteration. This element is added only if `time` flag is set to TRUE in 
#' config.}
#' 
#' \item{config}{The configuration used for processing.}
#' \item{call}{The function call itself.}
#' }
#' 
#' @examples 
#' 
#' fname <- system.file("extdata", "example_data.fa", 
#'                         package = "archR", mustWork = TRUE)
#' 
#' # Specifying 'dinuc' generates dinucleotide features
#' inputSeqsMat <- archR::prepare_data_from_FASTA(fasta_fname = fname,
#'     sinuc_or_dinuc = "dinuc")
#' 
#' inputSeqsRaw <- archR::prepare_data_from_FASTA(fasta_fname = fname, 
#'     raw_seq = TRUE)
#' 
#' # Set archR configuration
#' archRconfig <- archR::archR_set_config(
#'     parallelize = TRUE,
#'     n_cores = 2,
#'     n_iterations = 100,
#'     k_min = 1,
#'     k_max = 20,
#'     mod_sel_type = "stability",
#'     bound = 10^-8,
#'     inner_chunk_size = 100,
#'     flags = list(debug = FALSE, time = TRUE, verbose = TRUE,
#'         plot = FALSE)
#' )
#' 
#' # Run archR 
#' archRresult <- archR::archR(config = archRconfig,
#'                           seqs_ohe_mat = inputSeqsMat,
#'                           seqs_raw = inputSeqsRaw,
#'                           seqs_pos = seq(1,100,by=1),
#'                           threshold_itr = 2)
#' 
#'  
#' @export
archR <- function(config, seqs_ohe_mat, seqs_raw, seqs_pos = NULL,
                    threshold_itr = 3,
                    set_parsimony = c(FALSE, FALSE, FALSE),
                    set_ocollation = c(TRUE, FALSE, FALSE),
                    fresh = TRUE,
                    use_oc = NULL,
                    o_dir = NULL){
    ##
    flags <- config$flags
    plt <- config$flags$plotVerboseFlag
    dbg <- config$flags$debugFlag
    vrbs <- config$flags$verboseFlag
    tym <- config$flags$timeFlag
    chnksz <- config$innerChunkSize
    modSelType <- config$modSelType
    tol <- config$tol
    bound <- config$bound
    parallelize <- config$parallelize
    minSeqs <- config$minSeqs
    kFolds <- config$kFolds
    chkpnt <- config$checkpointing
    crs <- config$nCoresUse
    
    
    ## assert threshold_itr is a positive integer
    if(!threshold_itr > 0) {
        stop("Expecting threshold iteration to be numeric and > 0")
    }
    ## better, provide a summary function
    ## message(utils::str(config))
    archRStartTime <- Sys.time()
    ##
    if(!is.null(o_dir)){
        if(fresh){
            o_dir <- handle_dir_creation(o_dir, vrbs||dbg)
        }else{
            if(!dir.exists(o_dir)){
                stop(o_dir, " not found")
            }
        }
    }
    ##
    if(is.null(seqs_pos)){
        seqs_pos <- seq_len(Biostrings::width(seqs_raw[1]))
    }
    ##
    manage_o_dir(plt, o_dir) # this will stop if o_dir is NULL
    if(plt){
        plot_all_seqs_logo(seqs_raw, seqs_pos, dpath=o_dir)
    }
    ## Make checks for params in configuration
    .assert_archR_config(config, ncol(seqs_ohe_mat))
    .assert_archR_thresholdIteration(threshold_itr)
    ##
    ## ** To continue archR from an earlier run (further levels downstream)
    ## 1. Initializations of seqClustLabels, Factors etc should be
    ## appropriately handled.**
    ##
    ## Initialize sequence-cluster-labels and outerChunks
    seqsClustLabels <- rep("0", ncol(seqs_ohe_mat))
    seqsClustLabelsList <- vector("list", threshold_itr)
    clustFactors <- vector("list", threshold_itr)
    architectures <- vector("list", threshold_itr)
    if(tym) timeInfo <- vector("list", threshold_itr)
    ##
    test_itr <- 1
    ##----------------------------------------------------------------------
    if(fresh){
        outerChunksColl <- vector("list", 1)
        ## Set outerChunks for first iteration
        outerChunksColl[[1]] <- seq(ncol(seqs_ohe_mat))
        ##
    }else{
        .msg_pstr("Working on clusters from an earlier run", flg=vrbs)
        if(is.null(use_oc)){
            stop("'use_oc' should not be NULL")
        }else{
            .msg_pstr("OK", flg=dbg)
            ## use_oc is same as nxtOuterChunkColl
            seqsClustLabels <- .update_cluster_labels(seqsClustLabels,
                                                        use_oc)
            .assert_archR_OK_for_nextIteration(use_oc)
            outerChunksColl <- use_oc
            ##
        }
    }
    ##
    ##
    if(length(set_parsimony) < threshold_itr){
        set_parsimony <- rep(FALSE, threshold_itr)
        set_parsimony[length(set_parsimony)] <- TRUE
    }
    ##
    if(length(set_ocollation) < threshold_itr){
        set_ocollation <- rep(FALSE, threshold_itr)
        set_ocollation[1] <- TRUE
        set_ocollation[length(set_parsimony)] <- FALSE #Earlier was TRUE
    }
    
    #### Start cluster only once
    if(parallelize){
        cl <- parallel::makeCluster(crs, type = "FORK")
        parallel::setDefaultCluster(cl)
        .msg_pstr("Parallelization w/", crs, "cores", flg=dbg)
        .msg_pstr(cl, flg=dbg)
    }
    ####
    if(modSelType == "cv"){
        .msg_pstr("Model selection by cross-validation", flg=vrbs || dbg)
    }
    if(modSelType == "stability"){
        .msg_pstr("Model selection by factor stability", flg=vrbs || dbg)
        .msg_pstr("Bound:", bound, flg=vrbs || dbg)
    }
    

    ##
    ##-------------------------------------------------------------------------
    while (test_itr <= threshold_itr) {
        iterStartTime <- Sys.time()
        totOuterChunksColl <- length(outerChunksColl)
        .msg_pstr("=== Iteration", test_itr, 
            paste0("[", totOuterChunksColl, " chunk(s)]"), "===", flg=vrbs)
        nxtOuterChunksColl <- vector("list")
        seqsClustLabels <- rep("0", ncol(seqs_ohe_mat))
        intClustFactors <- NULL
        nClustEachOC <- rep(0, totOuterChunksColl)
        for (outerChunkIdx in seq_along(outerChunksColl)) {
            ##
            outerChunk <- outerChunksColl[[outerChunkIdx]]
            .msg_pstr(paste0("-----[Outer chunk ", outerChunkIdx, "/",
                        totOuterChunksColl, "]"), 
                paste0("[Size:", length(outerChunk), "]"), flg=vrbs)
            
            ## Make a decision to process based on size of chunk
            doNotProcess <- .decide_process_outer_chunk(minSeqs,
                                                        length(outerChunk),
                                                        kFolds)
            ##
            if (doNotProcess) {
                .msg_pstr("Decision: Skipping", flg=vrbs)
                ## TO-DO: Could write function to manipulate
                ## collatedClustAssignments
                ## NOTE: Control enters here only from second iteration onwards,
                ## otherwise statements here can fail. Because, clustFactors is
                ## assumed to be populated already -- in the first iteration,
                ## test_itr is 0, and indexing using test_itr would fail.
                collatedClustAssignments <- list(outerChunk)
                ## We access clustFactors that are set in the previous 
                ## iteration, when test_itr would be one less than its 
                ## current value.
                if (!is.null(intClustFactors)) {
                    intClustFactors <- 
                    cbind(intClustFactors, as.matrix(
                    clustFactors[[test_itr-1]]$basisVectors[, outerChunkIdx])
                    )
                } else {
                    intClustFactors <- as.matrix(
                    clustFactors[[test_itr-1]]$basisVectors[, outerChunkIdx]
                    )
                }
            } else {
                .msg_pstr("Decision: Processing", flg=vrbs)
                innerChunksColl <- .prepare_chunks(outerChunk,
                                            chnksz)
                ## Maintain these in a list, for collation later when
                ## all innerChunks in innerChunksColl have been processed
                globFactors <- vector("list", length(innerChunksColl))
                globClustAssignments <- vector("list", length(innerChunksColl))
                nClustEachIC <- rep(0, length(innerChunksColl))
                #################### INNER CHUNK FOR LOOP #####################
                for (innerChunkIdx in seq_along(innerChunksColl)) {
                    .msg_pstr(paste0("[Inner chunk ", innerChunkIdx, "/",
                        length(innerChunksColl), "]"), 
                        paste0("[Size:", 
                        length(innerChunksColl[[innerChunkIdx]]),"]"), flg=vrbs)
                    ##
                    ## Setting up sequences for the current chunk
                    this_seqsMat <-
                        seqs_ohe_mat[, innerChunksColl[[innerChunkIdx]]]
                    
                    if(test_itr == 1 ||
                        length(outerChunk) > 0.9*chnksz){
                        thisNMFResult <-
                            .handle_chunk_w_NMF2(innerChunkIdx,
                                        innerChunksColl,
                                        this_seqsMat,
                                        cgfglinear = TRUE,
                                        coarse_step = 10,
                                        monolinear = FALSE,
                                        askParsimony = set_parsimony[test_itr],
                                        doRegularize = FALSE,
                                        config, 
                                        o_dir, test_itr, 
                                        outerChunkIdx)
                    }else{
                        ##
                        thisNMFResult <-
                            .handle_chunk_w_NMF2(innerChunkIdx,
                                        innerChunksColl,
                                        this_seqsMat,
                                        cgfglinear = TRUE,
                                        coarse_step = 5,
                                        monolinear = TRUE,
                                        askParsimony = set_parsimony[test_itr],
                                        doRegularize = TRUE,
                                        config, 
                                        o_dir, test_itr, 
                                        outerChunkIdx)
                    }
                    .assert_archR_NMFresult(thisNMFResult)
                    globFactors[[innerChunkIdx]] <- 
                            thisNMFResult$forGlobFactors
                    globClustAssignments[[innerChunkIdx]] <-
                            thisNMFResult$forGlobClustAssignments
                    ##
                    nClustEachIC[innerChunkIdx] <- 
                            length(globClustAssignments[[innerChunkIdx]])
                } ## for loop over innerChunksColl ENDS here
                #################### INNER CHUNK FOR LOOP ######################
                ## We need globFactors, globClustAssignments.
                ## 
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
                    .collate_clusters2(globFactorsClustering,
                                        globClustAssignments)
                ## Collect number of clusters for each outer chunk
                nClustEachOC[outerChunkIdx] <- length(collatedClustAssignments)
            }  ## IfElse doNotProcess outer chunk ENDS
            ##
            .assert_archR_globClustAssignments(collatedClustAssignments)
            ## Assigning cluster labels can be done later, see below
            ## Collect (append) clusters at current level
            nxtOuterChunksColl <- append(nxtOuterChunksColl,
                                        collatedClustAssignments)
            ##
            chunksComplInfo <-
                paste0("-----[Outer chunk ", outerChunkIdx, "/", 
                            totOuterChunksColl, "] complete")
            iterComplInfo <- paste0("=====[Iteration ", test_itr, "] complete")
            currInfo <- paste0("current total factors: ",
                                ncol(intClustFactors))
            nextIterInfo <- paste0("Current total chunks for next iteration: ",
                                length(nxtOuterChunksColl))
            ##
            
            .msg_pstr(chunksComplInfo, flg=vrbs)
            if(outerChunkIdx == totOuterChunksColl) 
                .msg_pstr(iterComplInfo, flg=vrbs)
            .msg_pstr(chunksComplInfo, flg=dbg)
            .msg_pstr(currInfo, flg=dbg)
            .msg_pstr(nextIterInfo, flg=dbg)

        }  ## for loop over outerChunksCollection ENDS
        .msg_pstr("Managing clusters from outer chunk(s)", flg=vrbs)
        

        ############### Managing clusters from outer chunks
        intClustFactorsClustering <- NULL

        ## Can intClustFactors ever be NULL?
        
        .msg_pstr("totOuterChunks: ", totOuterChunksColl, flg=dbg)
        .msg_pstr(paste(nClustEachOC, collapse=","), flg=dbg)
        .msg_pstr(paste(nClustEachIC, collapse=","), flg=dbg)
        .msg_pstr("Collation this iter:", set_ocollation[test_itr],
                "; #InnerChunks: ", length(innerChunksColl), 
                "; #OuterChunks: ", totOuterChunksColl, flg=dbg)
        
        if (set_ocollation[test_itr] && (length(innerChunksColl) > 1 
            || totOuterChunksColl > 1)) {
            ## ^The second condition in the IF statement protects against
            ## performing HAC/collation when #inner chunks & the number of 
            ## outer chunks is 1.
            .msg_pstr("Decision for outer chunk collation: Yes", flg=dbg)
            ## Cluster the factors using hierarchical clustering
            ## Combinations considered:
            ## Cosine distance + Average Linkage
            ## Euclidean distance + Average Linkage
            ## modifiedNW distance + Average Linkage
            ## Euclidean distance + Complete Linkage
            ## Euclidean distance + Ward.D Linkage
            ## 
            ############## We use Complete linkage with Euclidean distance
            if(totOuterChunksColl > 1){
                .msg_pstr("meanClustersOC: ", ceiling(mean(nClustEachOC)), 
                    flg=dbg)
                ## For setting minClusters, note last iteration collated
                chkIdx <- seq_len(test_itr-1)
                if(any(set_ocollation[chkIdx])){
                    lastItrC <- tail(which(set_ocollation[chkIdx]), 1)
                    setMinClusters <- clustFactors[[lastItrC]]$nBasisVectors
                }else{
                    ## average clusters identified in each chunk of the 
                    ## first iteration
                    setMinClusters <- max(ceiling(mean(nClustEachOC[1])), 2) 
                }
            }else{
                ## When totOuterChunks is == 1, this is the first iteration
                ## Use the mean of nClustEachIC
                .msg_pstr("meanClustersIC: ", ceiling(mean(nClustEachIC)), 
                    flg=dbg)
                setMinClusters <- max(ceiling(mean(nClustEachIC)), 2)
            }
            ## regularize
            regIntClustFactors <- .regularizeMat(basisMat = intClustFactors,
                                                topN = 50)
            intClustFactorsClusteringEucCom <-
                .handle_clustering_of_factors(regIntClustFactors,
                                                clustMethod = "hc",
                                                linkage = "complete",
                                                distMethod = "cor",
                                                ## setting minClusters here
                                                ## can help to not undo the 
                                                ## clusters identified by NMF
                                                minClusters = setMinClusters,
                                                flags = flags)

            intClustFactors <- .get_factors_from_factor_clustering2(
                intClustFactorsClusteringEucCom, intClustFactors)
            ## Manage cluster assignments
            ## Call the hierarchical clustering version of collate_clusters
            ## This is named .collate_clusters2

            nxtOuterChunksColl <-
                .collate_clusters2(intClustFactorsClusteringEucCom,
                                    nxtOuterChunksColl, dbg)

            ## Updating cluster labels for sequences can be done later after 
            ## the clusters have been rearranged
            ## But we need this for rearrangement
            seqsClustLabels <- .update_cluster_labels(seqsClustLabels,
                                                        nxtOuterChunksColl)
            ##
        } else {
            .msg_pstr("Decision for outer chunk collation: No", flg=dbg)
            ## TO-DO: move this inside else block above?
            seqsClustLabels <- .update_cluster_labels(seqsClustLabels,
                                                        nxtOuterChunksColl)
        }
        ############### MANAGING CLUSTERS FROM OUTER CHUNK ENDS ################
        seqsClustLabelsList[[test_itr]] <- seqsClustLabels
        ##
        clustFactors[[test_itr]] <-
            .setup_clustFactors_for_archR_result(intClustFactors)
        ##
        .assert_archR_OK_for_nextIteration(nxtOuterChunksColl)
        outerChunksColl <- nxtOuterChunksColl
        if(plt){
            if(!is.null(o_dir)){
                intermediateResultsPlot(seqsClustLabelsList[[test_itr]],
                                    seqs_raw, pos_lab = seqs_pos,
                                    iter = test_itr, fname = o_dir,
                                    vrbs=vrbs||dbg)
            }
        }
        ##
        itrComplMsg <- paste("Iteration", test_itr, "completed")
        itrComplTime <- ""
        if(tym){
            complTime <- format(as.numeric(Sys.time() - iterStartTime,
                                    units = "mins"), digits = 3)
            itrComplTime <- paste(" in", complTime, "mins")
            timeInfo[[test_itr]] <- as.double.difftime(complTime)
        }
        .msg_pstr(itrComplMsg, itrComplTime, flg=vrbs)
        ##
        archRComplMsg <- paste("ellapsed since start")
        archRComplTime <- ""
        if(tym){
            complTime <- format(as.numeric(Sys.time() - archRStartTime,
                                        units = "mins"), digits = 3)
            archRComplTime <- paste(complTime, "mins ")
            .msg_pstr(archRComplTime, archRComplMsg, flg=vrbs)
        }
        ##
        ## test_itr was updated here
        ## write current iteration clustering to disk as RDS;
        ## good for checkpointing archR
        if(!is.null(o_dir) && chkpnt && test_itr != threshold_itr){
            .msg_pstr("Checkpointing set to TRUE. ",
                "Saving result of iteration ", test_itr, ".", flg=vrbs)
            curr_archRresult <- list(seqsClustLabels = seqsClustLabelsList,
                                        clustBasisVectors = clustFactors,
                                        rawSeqs = seqs_raw,
                                        config = config,
                                        call = match.call())
            rdsFilename <- paste0(o_dir, "archRresult_checkpoint",
                                    test_itr, ".rds")
            saveRDS(curr_archRresult, file=rdsFilename)
        }
        ##
        test_itr <- test_itr + 1
    } ## algorithm while loop ENDS
    ##
    if(tym){
        temp_res <- list(seqsClustLabels = seqsClustLabelsList,
                        clustBasisVectors = clustFactors,
                        rawSeqs = seqs_raw,
                        timeInfo = timeInfo,
                        config = config,
                        call = match.call())
    }else{
        temp_res <- list(seqsClustLabels = seqsClustLabelsList,
                        clustBasisVectors = clustFactors,
                        rawSeqs = seqs_raw,
                        config = config,
                        call = match.call())
    }
    ##
    decisionToCollate <- TRUE
    iterations <- length(clustFactors)
    if(iterations > 1){
        prevIterNFactors <- clustFactors[[iterations-1]]$nBasisVectors
        currIterNFactors <- clustFactors[[iterations]]$nBasisVectors
        if(currIterNFactors == prevIterNFactors){
            decisionToCollate <- FALSE
            .msg_pstr("Reordering decision: FALSE", flg=dbg)
        }
    }
    ## For setting minClusters, note last iteration collated
    if(any(set_ocollation)){
        lastItrC <- tail(which(set_ocollation), 1)
        setMinClustersFinal <- 
            temp_res$clustBasisVectors[[lastItrC]]$nBasisVectors
    }else{
        setMinClustersFinal <- 2 ## the default value
    }
    ##
    temp_res_reord <- collate_archR_result(temp_res,
                                    iter = threshold_itr,
                                    clust_method = "hc",
                                    aggl_method = "complete", ## or average?
                                    dist_method = "cor",
                                    minClusters = setMinClustersFinal,
                                    regularize = TRUE, 
                                    topn = 50,#floor(0.50*length(seqs_pos)),
                                    collate = decisionToCollate,
                                    flags = flags,
                                    enableSwitchSilToCH = FALSE)
    ## Print final stage output files to disk
    if(plt){
        if(!is.null(o_dir)){
            intermediateResultsPlot(seq_lab = temp_res_reord$seqsClustLabels,
                                    seqs_raw = seqs_raw, pos_lab = seqs_pos,
                                    iter = "Final", fname = o_dir,
                                    vrbs=vrbs||dbg)
        }
    }
    ##
    if(tym){
        temp_archRresult <- list(seqsClustLabels = seqsClustLabelsList,
                                clustBasisVectors = clustFactors,
                                clustSol = temp_res_reord,
                                rawSeqs = seqs_raw,
                                timeInfo = timeInfo,
                                config = config,
                                call = match.call())
    }else{
        temp_archRresult <- list(seqsClustLabels = seqsClustLabelsList,
                                clustBasisVectors = clustFactors,
                                clustSol = temp_res_reord,
                                rawSeqs = seqs_raw,
                                config = config,
                                call = match.call())
    }
    .assert_archRresult(temp_archRresult)
    ##
    ## Write result to disk as RDS file
    if(!is.null(o_dir)){
        rdsFilename <- paste0(o_dir, "archRresult.rds")
        saveRDS(temp_archRresult, file=rdsFilename)
        ##
        .msg_pstr("=== Result saved to location ===", flg=vrbs)
        .msg_pstr(rdsFilename, flg=vrbs)
    }
    ##
    ## At the moment, it is preferable to leave archR final iteration unordered.
    ## The final result clusters is computed with a particular setting (default)
    ## of collate_archR_result function. 
    ## 
    ## The user can choose the agglomeration method and the distMethod to 
    ## achieve suitable clustering results.
    ##
    ## Stop cluster
    .msg_pstr("Stopping cluster...",flg=dbg)
    if(parallelize) parallel::stopCluster(cl)
    ##
    archRComplMsg <- paste("archR exiting")
    archRComplTime <- ""
    if(tym){
        complTime <- format(as.numeric(Sys.time() - archRStartTime,
                                units = "mins"), digits = 3)
        archRComplTime <- paste(",", complTime, "mins")
    }
    .msg_pstr(archRComplMsg, archRComplTime, flg=vrbs)
    ##
    return(temp_archRresult)
} ## archR function ENDS
## =============================================================================
