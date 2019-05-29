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
#' @param seqs For representing the sequences as an image (TO-DO: May not be needed).
#' Default value is NULL when no image gets printed.
#' @param position_labels Labels of the positions in the sequences. Used for
#' visualization with function \code{viz_all_factors_in_combined_heatmaps_seqlogos}.
#'
#' @return nothing.
#'
#' @export
#' @importFrom stats quantile sd
plot_arch_for_clusters <- function(givenSamplesMatrix, givenFeaturesMatrix,
                                   nCluster, clustering_sol, seqs = NULL,
                                   position_labels = NA){
      # print(length(levels(as.factor(clustering_sol$clust_sol$cluster))))
      # print(levels(as.factor(clustering_sol$clust_sol$cluster)))
      if(!is.matrix(givenSamplesMatrix)){
            stop("givenSamplesMat not of type matrix")
      }
      if(sum(dim(givenSamplesMatrix)) == 2 && is.na(givenSamplesMatrix)){
            stop("Empty givenSamplesMat")
      }
      if(!is.matrix(givenFeaturesMatrix)){
            stop("givenFeaturesMat not of type matrix")
      }
      if(sum(dim(givenFeaturesMatrix)) == 2 && is.na(givenFeaturesMatrix)){
            stop("Empty givenFeaturesMat")
      }
      if(length(nCluster) > 1){
            stop("Expecting only one value for nCluster")
      }else if(nCluster < 2 ){
            stop("nCluster should be non-negative")
      }else if(length(levels(as.factor(clustering_sol$clust_sol$cluster))) != nCluster){
            stop("nCluster value and #Clusters in clustering_sol mismatch")
      }
      for (grp_ID in 1:nCluster){
            select_block <- givenSamplesMatrix[, clustering_sol$reordering_idx[[grp_ID]]]
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
            if(ncol(givenFeaturesMatrix) == nrow(mean_patt_in_block)){
                  feat_in_block <- givenFeaturesMatrix %*% mean_patt_in_block
            }else{
                  stop("Check dimensions of FeaturesMatrix and SamplesMatrix")
            }
            print(summary(feat_in_block))
            # sp_feat_in_block <- sparsify_mat(feat_in_block)
            # rescale
            feat_in_block <- (feat_in_block - mean(feat_in_block))/stats::sd(feat_in_block)
            print(summary(feat_in_block))
            print(as.numeric(quantile(feat_in_block, 0.95)))
            feat_in_block[feat_in_block < as.numeric(stats::quantile(feat_in_block, 0.95))] <- 10^-5
            # feat_in_block <- feat_in_block *  (norm(givenFeaturesMatrix, type="2")/norm(feat_in_block, type="2"))
            #
            print(summary(feat_in_block))
            # viz_all_factors_as_heatmap(feat_in_block)
            cluster_seqs <- seqs[ , clustering_sol$clust_sol$cluster == grp_ID]
            cat("#Sequences in Cluster ", grp_ID, ": ", clustering_sol$clust_sol$size[grp_ID], "\n")
            if(!is.null(seqs)){
                  represent_matrix_of_acgt( cluster_seqs,
                                            position_labels = position_labels,
                                            plot.title = paste0(clustering_sol$clust_sol$size[grp_ID],
                                                                " sequences in cluster ", grp_ID)
                                          )
            }
            #
            viz_all_factors_in_combined_heatmaps_seqlogos(feat_in_block, position_labels = position_labels, add_pseudo_counts = F)
      }


}
