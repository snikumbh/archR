#' @title Visualize the NMF basis vectors as a heatmap
#'
#' @description Each factor in the given features matrix is separately
#' represented as a heatmap. Since each feeature vector is a one-hot-encoded DNA
#'  sequence, it is reshaped to have 4 rows and then visualized.
#'
#' @inheritParams viz_basis_vectors_in_combined_heatmaps_seqlogos
#'
#' @return nothing.
#' @family visualization functions
#' @seealso \code{\link{viz_basis_vectors_as_seqlogo}} for plotting only as
#' sequence logos, \code{\link{viz_basis_vectors_in_combined_heatmaps_seqlogos}}
#' for plotting combined heatmaps and sequence logos.
#' @export
#'
viz_basis_vectors_as_heatmap <- function(feat_mat, pos_lab = NA,
                                        add_pseudo_counts = FALSE,
                                        pdf_name = NULL,
                                        sinuc_or_dinuc = "sinuc") {
    # Visualize all basis factors (expected as columns of the given features
    # matrix) as heatmaps
    if (!is.matrix(feat_mat)) {
        stop("feat_mat not of type matrix")
    }
    if (sum(dim(feat_mat)) == 2 && is.na(feat_mat)) {
        stop("Empty feat_mat")
    }

    if (sinuc_or_dinuc == "sinuc") {
        invisible(apply(feat_mat, MARGIN = 2, function(x) {
            pwm <- make_sinuc_PWMs(x, add_pseudo_counts = FALSE)
            p1 <- plot_ggheatmap(pwm_mat = pwm,
                                    pos_lab = pos_lab,
                                    pdf_name = pdf_name)
            base::suppressMessages(print(p1))
        }))
    } else if (sinuc_or_dinuc == "dinuc") {
        invisible(apply(feat_mat, MARGIN = 2, function(x) {
            pwm <- make_dinuc_PWMs(x, add_pseudo_counts = FALSE)
            p1 <- plot_ggheatmap(pwm_mat = pwm,
                                    pos_lab = pos_lab,
                                    pdf_name = pdf_name)
            base::suppressMessages(print(p1))
        }))
    }
}
