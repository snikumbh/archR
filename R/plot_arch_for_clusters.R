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
                                   nCluster, 
                                   clustering_sol, 
                                   seqs = NULL,
                                   position_labels = NA, 
                                   plotMethod = "custom", 
                                   add_pseudo_counts = F,
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
  } else if (nCluster < 2) {
      stop("nCluster should be non-negative")
  } #else if (length(levels(as.factor(clustering_sol$clust_sol$cluster))) != nCluster) {
    #  stop("nCluster value and #Clusters in clustering_sol mismatch")
  #}
  # cleaned-up code
  #
  dummy_samarth <- as.vector(rep(0,nrow(givenFeaturesMatrix)))
  architecture_features <- vector("list", ncol(givenFeaturesMatrix))
  architecture_features <- lapply(architecture_features, function(x){ x <- dummy_samarth})
  #
  for (grp_ID in 1:nCluster) {
    
      
      select_block <- givenSamplesMatrix[, clustering_sol$reordering_idx[[grp_ID]]]
      # Calculating means (since seqs along columns, we compute rowMeans)
      mean_patt_in_block <- as.matrix(rowMeans(select_block)) # asmatrix creates a column vector
      median_patt_in_block <- as.matrix(rowMedians(select_block))
      
      if (ncol(givenFeaturesMatrix) == nrow(mean_patt_in_block)) {
        meanFeat_in_block <- givenFeaturesMatrix %*% mean_patt_in_block
        medianFeat_in_block <- givenFeaturesMatrix %*% median_patt_in_block
        #
        
        #
        #
        new_meanFeat_in_block <- meanFeat_in_block
        new_medianFeat_in_block <- medianFeat_in_block
        new_meanFeat_in_block[meanFeat_in_block < as.numeric(stats::quantile(meanFeat_in_block, 0.75))] <- 10^-5
        new_medianFeat_in_block[medianFeat_in_block < as.numeric(stats::quantile(medianFeat_in_block, 0.75))] <- 10^-5
        #
        architecture_features[[grp_ID]] <- as.vector(t(meanFeat_in_block))
        #
      } else {
        stop("Check dimensions of FeaturesMatrix and SamplesMatrix")
      }
      
      # viz_all_factors_as_heatmap(meanFeat_in_block, position_labels = position_labels)
      if(clustering_sol$clustType == "kmeans"){
        cluster_seqs <- seqs[, clustering_sol$clust_sol$cluster == grp_ID]  
      } else if (clustering_sol$clustType == "hclust" || clustering_sol$clustType == "dbscan"){
        cluster_seqs <- seqs[, clustering_sol$clust_sol == grp_ID]
        
      }
      # else if(clustering_sol$clustType == "dbscan") {
      #   cluster_seqs <- seqs[, clustering_sol$clust_sol == grp_ID]
      #   
      # }
      
      #cat("#Sequences in Cluster ", grp_ID, ": ", clustering_sol$clust_sol$size[grp_ID], "\n")
      
      
      # if (!is.null(seqs)) {
      #   represent_matrix_of_acgt(cluster_seqs,
      #                            position_labels = position_labels,
      #                            plot.title = paste0(
      #                              clustering_sol$clust_sol$size[grp_ID],
      #                              " sequences in cluster "#, grp_ID
      #                            )
      #   )
      # }
      # Plotting mean features
      viz_all_factors_in_combined_heatmaps_seqlogos(meanFeat_in_block,
                                                    plotMethod = plotMethod,
                                                    position_labels = position_labels,
                                                    add_pseudo_counts = add_pseudo_counts,
                                                    sinuc_or_dinuc = sinuc_or_dinuc)
      # 
      # viz_all_factors_in_combined_heatmaps_seqlogos(meanFeat_in_block,
      #                                               plotMethod = "custom",
      #                                               position_labels = position_labels,
      #                                               add_pseudo_counts = add_pseudo_counts)
      
      
      
      # Plotting median features -- method=bits
      # viz_all_factors_in_combined_heatmaps_seqlogos(medianFeat_in_block,
      #                                               plotMethod = "bits",
      #                                               position_labels = position_labels,
      #                                               add_pseudo_counts = add_pseudo_counts)

      # viz_all_factors_in_combined_heatmaps_seqlogos(new_medianFeat_in_block,
      #                                               plotMethod = "bits",
      #                                               position_labels = position_labels,
      #                                               add_pseudo_counts = add_pseudo_counts)
      # # # Plotting median features -- method=custom
      # viz_all_factors_in_combined_heatmaps_seqlogos(medianFeat_in_block,
      #                                               plotMethod = "custom",
      #                                               position_labels = position_labels,
      #                                               add_pseudo_counts = add_pseudo_counts)
      # 
      # viz_all_factors_in_combined_heatmaps_seqlogos(new_medianFeat_in_block,
      #                                               plotMethod = "custom",
      #                                               position_labels = position_labels,
      #                                               add_pseudo_counts = add_pseudo_counts)
      
  }
  #
  # # Cross-correlating architectural features
  # print(architecture_features)
  # ccf_values <- matrix(rep(NA, nCluster*nCluster), byrow=T, nrow=nCluster)
  # for (grp_ID1 in 1:(nCluster-1)) {
  #   for (grp_ID2 in (grp_ID1+1):nCluster) {
  #     #print(class(architecture_features[[1]]))
  #     #ccf(architecture_features[[grp_ID1]],architecture_features[[grp_ID2]])
  #     #ccf_values_samarth <- ccf(architecture_features[[grp_ID1]], architecture_features[[grp_ID2]], type = "correlation", lag.max = 40)
  #     ccf_values_samarth <- forecast::Ccf(architecture_features[[grp_ID1]],
  #                                         architecture_features[[grp_ID2]],
  #                                         type = "correlation",
  #                                         lag.max = 5,
  #                                         plot = TRUE
  #                                         )
  #     #
  #     significance_threshold <- 2 * 1/(sqrt(length(architecture_features[[1]])))
  #     print(significance_threshold)
  #     #
  #     best_ccf <- max(ccf_values_samarth$acf)
  #     print(best_ccf)
  #     print(paste0(grp_ID1,"--",grp_ID2))
  #     print(stats::quantile(ccf_values_samarth$acf))
  #     print(stats::quantile(ccf_values_samarth$acf,0.9))
  #     if ( best_ccf < -1*significance_threshold || best_ccf > significance_threshold){
  # 
  #       #if()
  #       ccf_values[grp_ID1,grp_ID2] <- best_ccf
  #     }
  #   }
  # }
  # print("---")
  # #forecast::ggCcf(architecture_features, architecture_features, type = "correlation", lag.max = 20)
  # print(ccf_values)
  
  
  # # trials
  # for (grp_ID in 1:nCluster) {
  #   select_block <- givenSamplesMatrix[, clustering_sol$reordering_idx[[grp_ID]]]
  #   # Calculating means (since seqs along columns, we compute rowMeans)
  #   # print(select_block)
  #   mean_patt_in_block <- as.matrix(rowMeans(select_block)) # asmatrix creates a column vector
  #   # print(summary(mean_patt_in_block))
  #   # print(t(mean_patt_in_block))
  #   # print(sum(mean_patt_in_block))
  #   # divide by rowSums
  #   # print(summary(mean_patt_in_block/sum(mean_patt_in_block)))
  #   # mean_patt_in_block <- mean_patt_in_block/sum(mean_patt_in_block)
  #   # print(summary(mean_patt_in_block))
  #   # mean_patt_in_block <- (mean_patt_in_block - min(mean_patt_in_block)) /(max(mean_patt_in_block) - min(mean_patt_in_block))
  #   # print(sum(mean_patt_in_block))
  #   # print(summary( (mean_patt_in_block - min(mean_patt_in_block)) /(max(mean_patt_in_block) - min(mean_patt_in_block))   ))
  #   #
  #   if (ncol(givenFeaturesMatrix) == nrow(mean_patt_in_block)) {
  #     feat_in_block <- givenFeaturesMatrix %*% mean_patt_in_block
  #   } else {
  #     stop("Check dimensions of FeaturesMatrix and SamplesMatrix")
  #   }
  #   # print(summary(feat_in_block))
  #   # sp_feat_in_block <- sparsify_mat(feat_in_block)
  #   # rescale
  #   # feat_in_block <- (feat_in_block - mean(feat_in_block)) / stats::sd(feat_in_block)
  #   # print(summary(feat_in_block))
  #   # print(as.numeric(quantile(feat_in_block, 0.95)))
  #   # feat_in_block[feat_in_block < as.numeric(stats::quantile(feat_in_block, 0.50))] <- 10^-5
  #   # feat_in_block <- feat_in_block *  (norm(givenFeaturesMatrix, type="2")/norm(feat_in_block, type="2"))
  #   #
  #   # print(summary(feat_in_block))
  #   # viz_all_factors_as_heatmap(feat_in_block, position_labels = position_labels)
  #   cluster_seqs <- seqs[, clustering_sol$clust_sol$cluster == grp_ID]
  #   cat("#Sequences in Cluster ", grp_ID, ": ", clustering_sol$clust_sol$size[grp_ID], "\n")
  #   # #
  #   # decoded_seqs <- apply(cluster_seqs, 2, one_hot_decode)
  #   # decoded_seqs <- list(decoded_seqs)
  #   # #
  #   # p1 <- ggplot() + geom_logo(decoded_seqs, method = "bits", col_scheme = "nucleotide") +
  #   #   theme_logo() +
  #   #   scale_x_continuous(
  #   #     breaks = 1:(nrow(givenFeaturesMatrix)/4), labels = position_labels,
  #   #     expand = expand_scale(mult = c(0, 0))
  #   #   ) +
  #   #   theme(axis.text.x = element_text(size=rel(0.5), angle = 90, hjust = 1))
  #   # 
  #   #  # print(p1)
  #   
  #   if (!is.null(seqs)) {
  #     represent_matrix_of_acgt(cluster_seqs,
  #      position_labels = position_labels,
  #      plot.title = paste0(
  #        clustering_sol$clust_sol$size[grp_ID],
  #        " sequences in cluster "#, grp_ID
  #      )
  #     )
  #   }
  #   #
  #   # viz_all_factors_as_seqlogo(feat_in_block,
  #   #                             plotMethod = plotMethod,
  #   #                             position_labels = position_labels,
  #   #                             add_pseudo_counts = add_pseudo_counts)
  #   #
  #   viz_all_factors_in_combined_heatmaps_seqlogos(feat_in_block,
  #                                                 plotMethod = plotMethod,
  #                                                 position_labels = position_labels,
  #                                                 add_pseudo_counts = add_pseudo_counts)
  # }
  
  return (architecture_features)
}


get_seq_cluster_levels <- function(annClusters){
  print("fetching levels")
  clusterLevels <- levels(as.factor(annClusters))
  # interimBreaks <- vapply(1:length(clusterLevels), function(x){
  #   max(which(annClusters == clusterLevels[x]))
  # }, numeric(1))
  # seqClustBreaks <- c(1, interimBreaks)
  
  return(clusterLevels)
}

plot_arch_for_clusters_new <- function(tss.seqs_raw, archRresult, position_labels){
  
  
  # ggseqlogo::ggseqlogo(tss.seqs_raw[1:25])
  for(i in 1:2){#length(archRresult$clustFactors)){
    chosenLevelLabels <- collect_cluster_labels(archRresult$seqsClustLabels,
                                                choose_levels = i+1)
    chosenLevelLabels_sorted <- sort(chosenLevelLabels, index.return = TRUE)
    # print(chosenLevelLabels_sorted)
    fname <- paste0("samarth_trial_Level", i,"_architectures.pdf")
    cluster_levels <- get_seq_cluster_levels(chosenLevelLabels)
    # print(cluster_levels)
    print("print now? samarth")
    pdf(fname, width = 11, height = 1.5)
      samarth_invisible <- sapply(1:length(cluster_levels), function(x){
              relIdx <- which(chosenLevelLabels_sorted$x == cluster_levels[x])
              # print(relIdx)
              samarth_p <- ggseqlogo::ggseqlogo(tss.seqs_raw[relIdx])
              print(samarth_p)
              ##
              # samarth_p <- ggplot() + geom_logo(tss.seqs_raw[relIdx], 
              #                                   method = "custom", 
              #                                   seq_type = "dna") +
              #   theme_logo() +
              #   # #
              #   # ggplot2::scale_x_continuous(
              #   #   breaks = 1:length(position_labels), labels = position_labels,
              #   #   expand = expand_scale(mult = c(0, 0))
              #   # ) +
              #   # #
              #   # ggplot2::theme(axis.text.x = element_text(size = rel(0.5), angle = 90, hjust = 1),
              #   #                axis.text.y = element_text(size = rel(0.5))) 
              #
              # print(samarth_p)
              ##
                })
    dev.off()
  }
}
