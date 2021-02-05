#' @title Visualize the NMF basis vectors as a sequence logo
#'
#' @description Each factor in the given features matrix is separately
#' represented as a sequence logo. Since each feature vector is a
#' one-hot-encoded
#' DNA sequence, it is reshaped to have 4 rows and then visualized.
#'
#' @inheritParams viz_basis_vectors_in_combined_heatmaps_seqlogos
#'
#' @return nothing.
#' @family visualization functions
#' @export
#' @seealso \code{\link{viz_basis_vectors_as_heatmap}} for plotting only as
#' heatmap, \code{\link{viz_basis_vectors_in_combined_heatmaps_seqlogos}}
#' for plotting combined heatmaps and sequence logos.
#' @import ggplot2
#' @import ggseqlogo
#'
viz_basis_vectors_as_seqlogo <- function(feat_mat,
                                        method = "custom",
                                        pos_lab = NA,
                                        add_pseudo_counts = FALSE,
                                        pdf_name = NULL,
                                        sinuc_or_dinuc = "sinuc") {
    # Visualize all basis factors (expected as columns of the given features
    # matrix) as seqlogos
    if (!is.matrix(feat_mat)) {
        stop("feat_mat not of type matrix")
    }
    if (sum(dim(feat_mat)) == 2 && is.na(feat_mat)) {
        stop("Empty feat_mat")
    }
    invisible(apply(feat_mat, MARGIN = 2, function(x) {
        if (sinuc_or_dinuc == "dinuc") {
            dna_alphabet <- c("A", "C", "G", "T")
            dna_alphabet_dinuc <- do.call(paste0, expand.grid(dna_alphabet,
                                                                dna_alphabet))
            pwm <- make_dinuc_PWMs(as.matrix(x), add_pseudo_counts = FALSE)
        } else if (sinuc_or_dinuc == "sinuc") {
            pwm <- make_sinuc_PWMs(x, add_pseudo_counts = FALSE)
        }
        p1 <- plot_ggseqlogo(pwm_mat = pwm, method = method,
                                pos_lab = pos_lab,
                                pdf_name = pdf_name)
        base::suppressMessages(print(p1))
    }))
}
## =============================================================================
