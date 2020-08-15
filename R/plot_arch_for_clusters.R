#' @title Plot the sequence architectures of the sequence clusters.
#'
#' @description The different architectures characteristic of the different
#' sequence clusters are visualized. Therefore, the function takes as arguments
#' the samplesMatrix and featuresMatrix output by NMF, the number of clusters,
#' the clustering solution obtained and the position labels for the sequences.
#'
#'
#' @param givenSamplesMatrix The samples matrix resulting for NMF.
#' @param givenFeaturesMatrix The features matrix resulting from NMF.
#' @param nCluster The number of clusters.
#' @param clustering_sol The clustering solution as returned by
#' \code{get_clusters} function.
#' @param seqs For representing the sequences as an image
#' (TO-DO: May not be needed).
#' Default value is NULL when no image gets printed.
#' @param position_labels Labels of the positions in the sequences. Used for
#' visualization with function
#' \code{viz_basis_vectors_in_combined_heatmaps_seqlogos}.
#' @param add_pseudo_counts OK
#' @param sinuc_or_dinuc 'sinuc' or 'dinuc'
#' @param plotMethod 'custom' or 'bits', passed to ggseqlogo
#'
#' @return nothing.
#'
#' @export
#' @importFrom stats quantile sd
#' @importFrom BiocGenerics rowMeans
plot_arch_for_clusters <- function(givenSamplesMatrix, givenFeaturesMatrix,
                                    nCluster, clustering_sol, seqs = NULL,
                                    position_labels = NA, plotMethod = "custom",
                                    add_pseudo_counts = FALSE,
                                    sinuc_or_dinuc = "sinuc") {
    # print(length(levels(as.factor(clustering_sol$clust_sol$cluster))))
    # print(levels(as.factor(clustering_sol$clust_sol$cluster)))
    if (!is.matrix(givenSamplesMatrix)) {
        stop("givenSamplesMat not of type matrix")
    }
    if (sum(dim(givenSamplesMatrix)) == 2 && is.na(givenSamplesMatrix)) {
        stop("Empty givenSamplesMat")
    }
    if (!is.matrix(givenFeaturesMatrix)) {
        stop("givenFeaturesMat not of type matrix")
    }
    if (sum(dim(givenFeaturesMatrix)) == 2 && is.na(givenFeaturesMatrix)) {
        stop("Empty givenFeaturesMat")
    }
    if (length(nCluster) > 1) {
        stop("Expecting only one value for nCluster")
    } else if (nCluster < 1) {
        stop("nCluster should be non-negative")
    }
    # else if (length(levels(as.factor(clustering_sol$clust_sol$cluster)))
    # != nCluster) {
    # stop('nCluster value and #Clusters in clustering_sol mismatch') }
    # cleaned-up code
    dummy_samarth <- as.vector(rep(0, nrow(givenFeaturesMatrix)))
    architecture_features <- vector("list", ncol(givenFeaturesMatrix))
    architecture_features <- lapply(architecture_features, function(x) {
        x <- dummy_samarth
    })
    #
    for (grp_ID in seq_len(nCluster)) {


        select_block <-
            givenSamplesMatrix[, clustering_sol$reordering_idx[[grp_ID]]]
        # Calculating means (since seqs along columns, we compute rowMeans)
        mean_patt_in_block <-
            as.matrix(BiocGenerics::rowMeans(select_block))
        # asmatrix creates a column vector
        # median_patt_in_block <-
        #     as.matrix(matrixStats::rowMedians(select_block))

        if (ncol(givenFeaturesMatrix) == nrow(mean_patt_in_block)) {
            meanFeat_in_block <- givenFeaturesMatrix %*% mean_patt_in_block
            # medianFeat_in_block <- givenFeaturesMatrix %*% median_patt_in_block
            #

            #
            new_meanFeat_in_block <- meanFeat_in_block
            # new_medianFeat_in_block <- medianFeat_in_block
            new_meanFeat_in_block[meanFeat_in_block <
                                        as.numeric(stats::quantile(
                                            meanFeat_in_block,
                0.75))] <- 10^-5
            # new_medianFeat_in_block[medianFeat_in_block <
            #                             as.numeric(stats::quantile(
            #                                 medianFeat_in_block,
            #     0.75))] <- 10^-5
            #
            architecture_features[[grp_ID]] <- as.vector(t(meanFeat_in_block))
            #
        } else {
            stop("Check dimensions of FeaturesMatrix and SamplesMatrix")
        }

        # viz_basis_vectors_as_heatmap(meanFeat_in_block,
        # position_labels = position_labels)
        if (clustering_sol$clustType == "kmeans") {
            cluster_seqs <-
                seqs[, clustering_sol$clust_sol$cluster == grp_ID]
        } else if (clustering_sol$clustType == "hclust" ||
                    clustering_sol$clustType == "dbscan") {
            cluster_seqs <- seqs[, clustering_sol$clust_sol == grp_ID]

        }
        # else if(clustering_sol$clustType == 'dbscan') { cluster_seqs <- seqs[,
        # clustering_sol$clust_sol == grp_ID] }

        # cat('#Sequences in Cluster ', grp_ID, ': ',
        # clustering_sol$clust_sol$size[grp_ID], '\n')


        # if (!is.null(seqs)) { represent_matrix_of_acgt(cluster_seqs,
        # position_labels =
        # position_labels, plot.title =
        # paste0( clustering_sol$clust_sol$size[grp_ID], '
        # sequences in cluster '#, grp_ID ) ) } Plotting mean features
        viz_basis_vectors_in_combined_heatmaps_seqlogos(
            meanFeat_in_block,
            plotMethod = plotMethod,
            position_labels = position_labels,
            add_pseudo_counts = add_pseudo_counts,
            sinuc_or_dinuc = sinuc_or_dinuc)
        # viz_basis_vectors_in_combined_heatmaps_seqlogos(meanFeat_in_block,
        # plotMethod =
        # 'custom', position_labels = position_labels, add_pseudo_counts =
        # add_pseudo_counts)
    }

    return(architecture_features)
}


get_seq_cluster_levels <- function(annClusters) {
    print("fetching levels")
    clusterLevels <- levels(as.factor(annClusters))
    # interimBreaks <- vapply(1:length(clusterLevels), function(x){
    # max(which(annClusters == clusterLevels[x])) }, numeric(1))
    # seqClustBreaks <- c(1, interimBreaks)

    return(clusterLevels)
}

plot_arch_for_clusters_new <- function(tss.seqs_raw,
                                        list_of_elements,
                                        position_labels,
                                        xt_freq = 1,
                                        PDFfname = "archR_sequence_architectures.pdf") {
    if(!is.null(PDFfname)) {
        pdf(file=PDFfname, width = 11, height = 2)
    }
    # seqs_clusters_as_list <- get_seqs_clusters_in_a_list(archRresult_clust_labels)
    seqs_clusters_as_list <- list_of_elements
    cluster_lengths <- unlist(lapply(seqs_clusters_as_list, length))
    cumsums_of_cluster_lengths <- cumsum(cluster_lengths)
    cluster_names <- sort(as.character(seq_along(seqs_clusters_as_list)))
    for (i in seq_along(seqs_clusters_as_list)) {
        ##
        if(i > 1){
            startN <- 1 + cumsums_of_cluster_lengths[i-1]
        }else{
            startN <- 1
        }
        endN <- cumsums_of_cluster_lengths[i]
        plot_title <- paste0("(", i , "/", length(seqs_clusters_as_list), ") Arch `",
                             cluster_names[i], "': ",
                             length(seqs_clusters_as_list[[i]]),
                             " sequences (",  startN, "-",  endN, ")" )
        message(plot_title)
        ##
        samarth_p <- plot_ggseqlogo_of_seqs(seqs = tss.seqs_raw[ seqs_clusters_as_list[[i]] ],
                                            position_labels = position_labels,
                                            xt_freq = xt_freq,
                                            title = plot_title)
        ##
        suppressMessages(print(samarth_p))

    }
    if(!is.null(PDFfname)) {
        dev.off()
    }
}


#' @title Plot sequence logo of a collection of sequences
#'
#' @description A wrapper to ggseqlogo plotting. Given a collection of
#' sequences, this function plots the sequence logo.
#'
#' @param seqs Collection of sequences as a Biostrings::DNAStringSet object.
#' @param position_labels Numeric vector of labels for sequence positions.
#' This should be the samelength as the width of the given sequences.
#' @param xt_freq Specify the frequency of the x-axis ticks.
#' @param title The title for the plot
#'
#' @export
plot_ggseqlogo_of_seqs <- function(seqs, position_labels, xt_freq = 1,
                                       title = "Title"){

    nPos <- length(position_labels)
    xtick_cal <- seq(0, nPos, by = xt_freq)
    xtick_cal[1] <- 1
    xtick_cal[length(xtick_cal)] <- nPos

    samarth_p <-
        ggseqlogo::ggseqlogo(
            as.character(seqs),
            seq_type = "dna",
            method = "bits"
        ) +
        ggplot2::theme_linedraw() +
        ggplot2::theme(axis.text.x = element_text(size = rel(0.7),
                                                  angle = 90,
                                                  hjust = 1),
                       axis.text.y = element_text(size = rel(0.8)),
                       panel.grid = element_blank()
        ) +
        # ggplot2::scale_x_continuous(breaks = seq_len(length(position_labels)),
        #                             labels = position_labels,
        #                             expand = expansion(mult = c(0, 0))) +
        ## Add additional bold tick labels
        ggplot2::scale_x_continuous(breaks = xtick_cal,
                                    labels = position_labels[xtick_cal],
                                    expand = expansion(mult = c(0, 0))) +

        ggplot2::ylim(0.0, 2.0) +
        ggplot2::ggtitle(title)
        message("Plot title:", title)

    samarth_p
}
