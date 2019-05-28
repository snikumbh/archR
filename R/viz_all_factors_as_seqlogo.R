#' @title Visualize the factors in features matrix as a sequence logo.
#'
#' @description Each factor in the given features matrix is separately
#' represented as a sequence logo. Since each feature vector is a one-hot-encoded
#' DNA sequence, it is reshaped to have 4 rows and then visualized.
#'
#' @inheritParams viz_all_factors_in_combined_heatmaps_seqlogos
#'
#' @return nothing.
#'
#' @export
#'
#' @import ggplot2
#' @import ggseqlogo
#'
#' @examples
#'
#'
viz_all_factors_as_seqlogo <- function(featuresMatrix, position_labels=NA, add_pseudo_counts = F, savePDFfilename=NULL){
      # Visualize all basis factors (expected as columns of the given features matrix)
      # as seqlogos
      if(!is.matrix(featuresMatrix)){
        stop("featuresMatrix not of type matrix")
      }
      if(sum(dim(featuresMatrix)) == 2 && is.na(featuresMatrix)){
        stop("Empty featuresMatrix")
      }
      invisible(apply(featuresMatrix, MARGIN = 2, function(x){
        pwm <- make_sinuc_PWMs(x, add_pseudo_counts = F)
        p1 <- plot_ggseqlogo(pwmMat=pwm,
                             position_labels=position_labels,
                             savePDFfilename=savePDFfilename
        )
        print(p1)
      }
      ))
}
