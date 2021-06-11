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
#' @param total_itr Numeric. Specify the number of iterations to perform. This
#' should be greater than zero. Default is NULL.
#' @param set_parsimony Logical vector. A logical vector of length `total_itr`
#' specifying if model selection by cross-validation should prefer parsimonious
#' solutions. Not required when stability-based model selection is chosen. TRUE
#' denotes parsimony is followed, FALSE otherwise.
#' @param set_ocollation Logical vector. A logical vector of length `total_itr`
#' specifying for every iteration of archR if collation of clusters from
#' outer chunks should be performed. TRUE denotes clusters are collated,
#' FALSE otherwise.
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
#' \code{\link{get_seqs_clust_list}} function), or simply clusters from any
#' other method.
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
#'                           total_itr = 2,
#'                           set_ocollation = c(TRUE, FALSE))
#'
#'
#' @export
archR <- function(config, seqs_ohe_mat, seqs_raw, seqs_pos = NULL,
                    total_itr = NULL,
                    set_parsimony = NULL,
                    set_ocollation = NULL,
                    fresh = TRUE,
                    use_oc = NULL,
                    o_dir = NULL){
    ##
    archRStartTime <- Sys.time()
    ##
    flags <- config$flags
    plt <- config$flags$plotVerboseFlag
    dbg <- config$flags$debugFlag
    vrbs <- config$flags$verboseFlag
    tym <- config$flags$timeFlag
    chnksz <- config$innerChunkSize
    modSelType <- config$modSelType
    bound <- config$bound
    parallelize <- config$parallelize
    minSeqs <- config$minSeqs
    kFolds <- config$kFolds
    chkpnt <- config$checkpointing
    crs <- config$nCoresUse

    cli::cli_rule(left="Setting up")

    setup_ans <- perform_setup(config, total_itr, o_dir, fresh,
                    seqs_pos, seqs_raw, seqs_ohe_mat, set_parsimony,
                    set_ocollation)
    seqs_pos <- setup_ans$seqs_pos
    o_dir <- setup_ans$o_dir



    ##
    ## ** To continue archR from an earlier run (further levels downstream)
    ## 1. Initializations of seqClustLabels, Factors etc should be
    ## appropriately handled.**
    ##
    ## Initialize sequence-cluster-labels and outerChunks
    seqsClustLabels <- rep("0", ncol(seqs_ohe_mat))
    seqsClustLabelsList <- vector("list", total_itr)
    clustFactors <- vector("list", total_itr)
    architectures <- vector("list", total_itr)
    timeInfo <- lapply(seq_len(total_itr), function(x){NA})
    ##
    test_itr <- 1
    ##--------------------------------------------------------------------------
    ## Set outerChunks for first iteration in fresh or non-fresh case
    if(fresh){
        outerChunksColl <- vector("list", 1)
        outerChunksColl[[1]] <- seq(ncol(seqs_ohe_mat))
    }else{
        cli::cli_alert_info("Working on clusters from an earlier run")
        if(is.null(use_oc)) stop("'use_oc' should not be NULL")
        ## use_oc is same as nxtOuterChunkColl
        seqsClustLabels <- .update_cluster_labels(seqsClustLabels, use_oc)
        .assert_archR_OK_for_nextIteration(use_oc)
        outerChunksColl <- use_oc
    }
    ##

    ##-------------------------------------------------------------------------
    while (test_itr <= total_itr) {
        iterStartTime <- Sys.time()
        totOuterChunksColl <- length(outerChunksColl)
        ##
        cli::cli_h1(c("Iteration {test_itr} of {total_itr} ",
            "[{totOuterChunksColl} chunk{?s}]"))
        ##
        nxtOuterChunksColl <- vector("list")
        seqsClustLabels <- rep("0", ncol(seqs_ohe_mat))
        intClustFactors <- NULL
        nClustEachOC <- rep(0, totOuterChunksColl)
        for (outerChunkIdx in seq_along(outerChunksColl)) {
            ##
            outerChunk <- outerChunksColl[[outerChunkIdx]]
            lenOC <- length(outerChunk)
            ##
            cli::cli_h2(c("Outer chunk {outerChunkIdx} ",
                "of {totOuterChunksColl} [Size: {lenOC}]"))
            ## Make a decision to process based on size of chunk
            doNotProcess <- .decide_process_outer_chunk(minSeqs, lenOC, kFolds)
            ##
            if (doNotProcess) {
                ##
                cli::cli_alert_info("Decision: Skipping")
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
                intClustFactors <-
                    cbind(intClustFactors, as.matrix(
                    clustFactors[[test_itr-1]]$basisVectors[, outerChunkIdx]))
                ##
            } else {
                innerChunksColl <- .prepare_chunks(outerChunk, chnksz)
                ## Maintain these in a list, for collation later when
                ## all innerChunks in innerChunksColl have been processed
                icResult <- process_innerChunk(test_itr, innerChunksColl,
                                    config, lenOC, seqs_ohe_mat, set_parsimony,
                                    o_dir, outerChunkIdx)
                globFactors <- icResult$globFactors
                globClustAssignments <- icResult$globClustAssignments
                nClustEachIC <- icResult$nClustEachIC

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
                tempClustFactors <- .get_factors_from_factor_clustering2(
                            globFactorsClustering, globFactorsMat)
                intClustFactors <- cbind(intClustFactors, tempClustFactors)
                ##
                ## Manage collated cluster assignments
                collatedClustAssignments <- .collate_clusters2(
                                globFactorsClustering, globClustAssignments)
                ## Collect number of clusters for each outer chunk
                nClustEachOC[outerChunkIdx] <- length(collatedClustAssignments)
            }  ## IfElse doNotProcess outer chunk ENDS
            ##
            .assert_archR_globClustAssignments(collatedClustAssignments)
            ## Assigning cluster labels can be done later, see below
            ## Collect (append) clusters at current level
            nxtOuterChunksColl <-
                append(nxtOuterChunksColl, collatedClustAssignments)
            ##
            cli::cli_alert_success(c("{outerChunkIdx} of {totOuterChunksColl} ",
                "outer chunk{?s} complete"))

            ##
            currInfo <- paste0("current total factors: ",
                                ncol(intClustFactors))
            nextIterInfo <- paste0("Current total chunks for next iteration: ",
                                length(nxtOuterChunksColl))
            .msg_pstr(c(currInfo, "\n", nextIterInfo), flg=dbg)
            ##
            if(outerChunkIdx == totOuterChunksColl) {
                cli::cli_alert_success(c("{test_itr} of ",
                    "{total_itr} iteration{?s} complete"))
            }
        }  ## for loop over outerChunksCollection ENDS
        .msg_pstr("Managing clusters from outer chunk(s)", flg=dbg)

        ############### Managing clusters from outer chunks
        intClustFactorsClustering <- NULL

        ## Can intClustFactors ever be NULL?
        #### Debugging print messages ####
        .msg_pstr("totOuterChunks: ", totOuterChunksColl, flg=dbg)
        .msg_pstr(paste(nClustEachOC, collapse=","), flg=dbg)
        .msg_pstr(paste(nClustEachIC, collapse=","), flg=dbg)
        .msg_pstr("Collation this iter:", set_ocollation[test_itr],
                "; #InnerChunks: ", length(innerChunksColl),
                "; #OuterChunks: ", totOuterChunksColl, flg=dbg)
        ####
        if (set_ocollation[test_itr] && (length(innerChunksColl) > 1
            || totOuterChunksColl > 1)) {
            ## ^The second condition in the IF statement protects against
            ## performing HAC/collation when #inner chunks & the number of
            ## outer chunks is 1.
            .msg_pstr("Decision for outer chunk collation: Yes", flg=dbg)
            ## Cluster the factors using hierarchical clustering
            setMinClusters <- keepMinClusters(set_ocollation, temp_res,
                                totOuterChunksColl, test_itr,
                                nClustEachOC, nClustEachIC, dbg,
                                clustFactors, stage = NULL)

            ## regularize basis matrix
            regIntClustFactors <- .regularizeMat(basisMat = intClustFactors,
                                                topN = 50)
            ## setting minClusters here can help to not undo the clusters
            ## identified by NMF
            intClustFactorsClusteringEucCom <-
                .handle_clustering_of_factors(regIntClustFactors,
                                                clustMethod = "hc",
                                                linkage = "complete",
                                                distMethod = "cor",
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
            ## the clusters have been rearranged. But we need this for
            ## rearrangement
            seqsClustLabels <- .update_cluster_labels(seqsClustLabels,
                                                        nxtOuterChunksColl)
            ##
        } else {
            .msg_pstr("Decision for outer chunk collation: No", flg=dbg)
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

        ## write current iteration clustering to disk as RDS;
        ## good for checkpointing archR
        if(chkpnt) save_checkpoint(o_dir, test_itr, total_itr,
                                seqsClustLabelsList, clustFactors,
                                seqs_raw, config, call = match.call())
        ##
        if(tym) {
            timeInfo[[test_itr]] <- show_ellapsed_time(
                use_str =  paste("Iteration", test_itr, "completed: "),
                use_time = iterStartTime)
            foo <- show_ellapsed_time(use_time = archRStartTime)
        }
        test_itr <- test_itr + 1
        ##
    } ## algorithm while loop ENDS
    ##
    temp_res <- list(seqsClustLabels = seqsClustLabelsList,
                    clustBasisVectors = clustFactors,
                    rawSeqs = seqs_raw,
                    timeInfo = timeInfo,
                    config = config,
                    call = match.call())
    ##

    decisionToCollate <- decisionToCollate(clustFactors, dbg=dbg)
    setMinClustersFinal <- keepMinClusters(set_ocollation, temp_res,
                                           totOuterChunksColl =
                                               totOuterChunksColl, dbg = dbg,
                                           nClustEachIC = nClustEachIC,
                                           test_itr = test_itr -1,
                                           stage="Final")
    ##
    temp_res_reord <- collate_archR_result(temp_res,
                                    iter = total_itr,
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
    if(plt && !is.null(o_dir)){
        intermediateResultsPlot(seq_lab = temp_res_reord$seqsClustLabels,
                                seqs_raw = seqs_raw, pos_lab = seqs_pos,
                                iter = "Final", fname = o_dir,
                                vrbs=vrbs||dbg)
    }
    ##
    temp_archRresult <- list(seqsClustLabels = seqsClustLabelsList,
                            clustBasisVectors = clustFactors,
                            clustSol = temp_res_reord,
                            rawSeqs = seqs_raw,
                            timeInfo = timeInfo,
                            config = config,
                            call = match.call())
    ##
    .assert_archRresult(temp_archRresult)
    ## Write result to disk as RDS file
    save_final_result(o_dir, temp_archRresult)

    ##
    ## It is preferable to leave archR final iteration unordered.
    ## The final result clusters is computed with a particular setting (default)
    ## of collate_archR_result function.
    ##
    ## The user can choose the agglomeration method and the distMethod to
    ## achieve suitable clustering results.

    ## Stop cluster
    if(parallelize) parallel::stopCluster(setup_ans$cl)
    ##
    if(tym){
        complTime1 <- Sys.time() - archRStartTime
        cli::cli_rule(c("archR exiting {prettyunits::pretty_dt(complTime1)}"))
    }
    ##
    return(temp_archRresult)
} ## archR function ENDS
## =============================================================================
