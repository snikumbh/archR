#' @title Visualize the NMF basis vectors in a paired heatmap and sequence logo
#' plot
#'
#' @description The given features matrix (with 4 rows) is represented as a
#' heatmap followed by a sequence logo where the positions are aligned for
#' better visualization.
#'
#' @param featuresMatrix The features matrix output from NMF. Expected
#' dimensionality: number of columns represent the number of factors from NMF,
#' and the number of rows is 4 times the length of the sequences in the
#' collection.
#' @param plotMethod For \code{ggseqlogo}; either of "custom", "bits", or
#' "probability".
#' Default is "custom".
#' @param position_labels Labels of the positions in the sequences.
#' @param add_pseudo_counts Boolean, taking values TRUE/T or FALSE/F, default
#' set
#'  to FALSE. Setting it to TRUE will enable adding pseudo-counts to the
#'  features
#'  matrix.
#' @param savePDFfilename Name of the file which will be saved as PDF
#' (also provide the extension).
#' @param sinuc_or_dinuc "sinuc" or "dinuc"
#'
#' @return nothing
#'
#' @export
#' @family visualization functions
#' @seealso \code{\link{viz_basis_vectors_as_heatmap}} for plotting only as
#' heatmaps, \code{\link{viz_basis_vectors_as_seqlogo}} for plotting only as
#' sequence logos.
#' @import ggplot2
#' @import ggseqlogo
viz_basis_vectors_in_combined_heatmaps_seqlogos <-
    function(featuresMatrix, plotMethod = "custom", position_labels = NA,
                add_pseudo_counts = FALSE, savePDFfilename = NULL,
                sinuc_or_dinuc = "sinuc") {
        ##
    if (!is.matrix(featuresMatrix)) {
        stop("featuresMatrix not of type matrix")
    }
    if (sum(dim(featuresMatrix)) == 2 && is.na(featuresMatrix)) {
        stop("Empty featuresMatrix")
    }
    invisible(apply(featuresMatrix, MARGIN = 2, function(x) {
        if (sinuc_or_dinuc == "dinuc") {
            dna_alphabet <- c("A", "C", "G", "T")
            dna_alphabet_dinuc <- do.call(paste0, expand.grid(dna_alphabet,
                                                                dna_alphabet))
            pwm <- make_dinuc_PWMs(as.matrix(x), add_pseudo_counts = FALSE,
                                scale = FALSE)
        } else if (sinuc_or_dinuc == "sinuc") {
            pwm <- make_sinuc_PWMs(x, add_pseudo_counts = FALSE, scale = FALSE)
        }
        ## Heatmap on top
        p1 <- plot_ggheatmap(
            pwmMat = pwm,
            position_labels = position_labels,
            savePDFfilename = savePDFfilename
        )
        p1 + theme(
            plot.margin = grid::unit(c(0, 0, 0, 0), "mm")
        )
        ## Seqlogo below
        p2 <- plot_ggseqlogo(
            pwmMat = pwm,
            plotMethod = plotMethod,
            position_labels = position_labels,
            savePDFfilename = savePDFfilename
        )
        ## Make adjustments for alignment
        p2 + theme(
            plot.margin = grid::unit(c(0, 0, 0, 0), "mm")
        )
        final_p <- ggpubr::ggarrange(p1, p2, ncol = 1, heights = c(1,1),
                            widths = c(1,1), align = 'v')
        ##
        if (!is.null(savePDFfilename)) {
            if (file.exists(savePDFfilename)) {
                warning("File exists, will overwrite")
            }
            ggplot2::ggsave(final_p, device = "pdf", width = 20, height = 2.5)
        } else {
            base::suppressMessages(print(final_p))
        }
    }))
}
## =============================================================================
