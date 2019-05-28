#' @title Plot the sequence architectures of the sequence clusters.
#'
#' @description The different architectures characteristic of the different
#' sequence clusters are visualized. Therefore, the function takes as arguments
#' the samplesMatrix and featuresMatrix output by NMF, the number of clusters,
#' the clustering solution obtained and the position labels for the sequences.
#'
#'
#' @param givenSamplesMat The samples matrix resulting for NMF.
#' @param givenFeaturesMatrix The features matrix resulting from NMF.
#' @param nCluster The number of clusters.
#' @param clustering_sol The clustering solution as returned by
#' \code{get_clusters} function.
#' @param seqs For representing the sequences as an image (TO-DO: May not be needed).
#' @param position_labels Labels of the positions in the sequences. Used for
#' visualization with function \code{viz_all_factors_in_combined_heatmaps_seqlogos}.
#'
#' @return nothing.
#'
#' @export
#' @examples
#'
plot_arch_for_clusters <- function(givenSamplesMat, givenFeaturesMatrix,
                                   nCluster, clustering_sol, seqs,
                                   position_labels){

  for (grp_ID in 1:nCluster){
        select_block <- givenSamplesMat[, clustering_sol$reordering_idx[[grp_ID]]]
        # Calculating means (since seqs along columns, we compute rowMeans)
        # print(select_block)
        mean_patt_in_block <- as.matrix(rowMeans(select_block)) # asmatrix creates a column vector
        # print(summary(mean_patt_in_block))
        # print(t(mean_patt_in_block))
        # print(sum(mean_patt_in_block))
        # divide by rowSums
        # print(summary(mean_patt_in_block/sum(mean_patt_in_block)))
        # mean_patt_in_block <- mean_patt_in_block/sum(mean_patt_in_block)
        # print(summary(mean_patt_in_block))
        #mean_patt_in_block <- (mean_patt_in_block - min(mean_patt_in_block)) /(max(mean_patt_in_block) - min(mean_patt_in_block))
        # print(sum(mean_patt_in_block))
        # print(summary( (mean_patt_in_block - min(mean_patt_in_block)) /(max(mean_patt_in_block) - min(mean_patt_in_block))   ))
        #
        feat_in_block <- givenFeaturesMatrix %*% mean_patt_in_block
        print(summary(feat_in_block))
        # sp_feat_in_block <- sparsify_mat(feat_in_block)
        # rescale
        feat_in_block <- (feat_in_block - mean(feat_in_block))/sd(feat_in_block)
        print(summary(feat_in_block))
        print(as.numeric(quantile(feat_in_block, 0.95)))
        feat_in_block[feat_in_block < as.numeric(quantile(feat_in_block, 0.95))] <- 10^-5
        # feat_in_block <- feat_in_block *  (norm(givenFeaturesMatrix, type="2")/norm(feat_in_block, type="2"))
        #
        print(summary(feat_in_block))
        # viz_all_factors_as_heatmap(feat_in_block)
        cluster_seqs <- seqs[ , clustering_sol$clust_sol$cluster == grp_ID]
        cat("#Sequences in Cluster ", grp_ID, ": ", clustering_sol$clust_sol$size[grp_ID], "\n")
        # represent_matrix_of_acgt( cluster_seqs,
        #                           positions,
        #                           plot.title = paste0(clustering_sol$clust_sol$size[grp_ID], " sequences in cluster ", grp_ID)
        #                         )
        #
        viz_all_factors_in_combined_heatmaps_seqlogos(feat_in_block, position_labels = positions_labels, add_pseudo_counts = F)
  }


}
