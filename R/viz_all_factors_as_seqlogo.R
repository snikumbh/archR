#' @title Visualize the factors in features matrix as a sequence logo
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
viz_all_factors_as_seqlogo <- function(featuresMatrix, plotMethod = "custom", position_labels = NA, 
    add_pseudo_counts = F, savePDFfilename = NULL, sinuc_or_dinuc = "sinuc") {
    # Visualize all basis factors (expected as columns of the given features matrix) as
    # seqlogos
    if (!is.matrix(featuresMatrix)) {
        stop("featuresMatrix not of type matrix")
    }
    if (sum(dim(featuresMatrix)) == 2 && is.na(featuresMatrix)) {
        stop("Empty featuresMatrix")
    }
    invisible(apply(featuresMatrix, MARGIN = 2, function(x) {
        if (sinuc_or_dinuc == "dinuc") {
            dna_alphabet <- c("A", "C", "G", "T")
            dna_alphabet_dinuc <- do.call(paste0, expand.grid(dna_alphabet, dna_alphabet))
            pwm <- make_dinuc_PWMs(as.matrix(x), add_pseudo_counts = F)
        } else if (sinuc_or_dinuc == "sinuc") {
            pwm <- make_sinuc_PWMs(x, add_pseudo_counts = F)
        }
        
        p1 <- plot_ggseqlogo(pwmMat = pwm, plotMethod = plotMethod, position_labels = position_labels, 
            savePDFfilename = savePDFfilename)
        base::suppressMessages(print(p1))
    }))
}
