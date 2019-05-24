viz_all_factors_as_heatmap <- function(featuresMatrix, position_labels=NA, add_pseudo_counts = F, savePDFfilename=NULL){
      # Visualize all basis factors (expected as columns of the given features matrix)
      # as heatmaps
      #
      invisible(apply(featuresMatrix, MARGIN = 2, function(x){
        pwm <- make_sinuc_PWMs(x, add_pseudo_counts = F)
        p1 <- plot_ggheatmap(pwmMat=pwm, 
                             position_labels=positions, 
                             savePDFfilename=savePDFfilename
        )
        print(p1)
      }
      ))
}