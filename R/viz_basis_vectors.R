#' @title Visualize the NMF basis vectors in a paired heatmap and sequence logo
#' plot
#'
#' @description The given features matrix (with 4 rows) is represented as a
#' heatmap followed by a sequence logo where the positions are aligned for
#' better visualization.
#'
#' @param feat_mat The features matrix output from NMF. Expected
#' dimensionality: number of columns represent the number of factors from NMF,
#' and the number of rows is 4 times the length of the sequences in the
#' collection.
#' @param method For \code{ggseqlogo}; either of "custom", "bits", or
#' "probability".
#' Default is "custom".
#' @param pos_lab Labels of the positions in the sequences.
#' @param add_pseudo_counts Boolean, taking values TRUE/T or FALSE/F, default
#' set
#'  to FALSE. Setting it to TRUE will enable adding pseudo-counts to the
#'  features
#'  matrix.
#' @param pdf_name Name of the file which will be saved as PDF
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
    function(feat_mat, method = "custom", pos_lab = NA, 
        add_pseudo_counts = FALSE, pdf_name = NULL,
                sinuc_or_dinuc = "sinuc") {
    ##
    check_vars2(feat_mat)
    # if(!requireNamespace("ggpubr", quietly = TRUE)){
    #     stop("Please install the R package 'ggpubr' to use this function")
    # }else{
    #     if (!is.matrix(feat_mat)) {
    #         stop("feat_mat not of type matrix")
    #     }
    #     if (sum(dim(feat_mat)) == 2 && is.na(feat_mat)) {
    #         stop("Empty feat_mat")
    #     }
    invisible(apply(feat_mat, MARGIN = 2, function(x) {
    if (sinuc_or_dinuc == "dinuc") {
        dna_alphabet <- c("A", "C", "G", "T")
        dna_alphabet_dinuc <- do.call(paste0, 
                        expand.grid(dna_alphabet, dna_alphabet))
        pwm <- make_dinuc_PWMs(as.matrix(x), add_pseudo_counts = FALSE,
                            scale = FALSE)
    } else if (sinuc_or_dinuc == "sinuc") {
        pwm <- make_sinuc_PWMs(x, add_pseudo_counts = FALSE, 
                            scale = FALSE)
    }
    ## Heatmap on top
    p1 <- plot_ggheatmap(pwm_mat = pwm, pos_lab = pos_lab,
                    pdf_name = pdf_name)
    p1 + theme(
        plot.margin = grid::unit(c(0, 0, 0, 0), "mm")
    )
    ## Seqlogo below
    p2 <- plot_ggseqlogo(pwm_mat = pwm, method = method,
        pos_lab = pos_lab, pdf_name = pdf_name
    )
    ## Make adjustments for alignment
    p2 + theme(plot.margin = grid::unit(c(0, 0, 0, 0), "mm"))
    final_p <- ggpubr::ggarrange(p1, p2, ncol = 1, heights = c(1,1),
                        widths = c(1,1), align = 'v')
    ##
    if (!is.null(pdf_name)) {
        if (file.exists(pdf_name)) {
            warning("File exists, will overwrite", immediate. = TRUE)
        }
        ggplot2::ggsave(final_p, device = "pdf", width = 20, 
                            height = 2.5)
    } else {
        base::suppressMessages(print(final_p))
    }
}))
    
}
## =============================================================================


check_vars2 <- function(feat_mat){
    if(!requireNamespace("ggpubr", quietly = TRUE)){
        stop("Please install the R package 'ggpubr' to use this function")
    }else{
        if (!is.matrix(feat_mat)) {
            stop("feat_mat not of type matrix")
        }
        if (sum(dim(feat_mat)) == 2 && is.na(feat_mat)) {
            stop("Empty feat_mat")
        }
    }
}
