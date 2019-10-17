#' @title Visualize the factors in features matrix as a heatmap plot
#'
#' @description Each factor in the given features matrix is separately
#' represented as a heatmap. Since each feeature vector is a one-hot-encoded DNA
#'  sequence, it is reshaped to have 4 rows and then visualized.
#'
#' @inheritParams viz_basis_vectors_in_combined_heatmaps_seqlogos
#'
#' @return nothing.

#' @export
#'
viz_basis_vectors_as_heatmap <- function(featuresMatrix, position_labels = NA,
                                        add_pseudo_counts = FALSE,
                                        savePDFfilename = NULL,
                                        sinuc_or_dinuc = "sinuc") {
    # Visualize all basis factors (expected as columns of the given features
    # matrix) as heatmaps
    if (!is.matrix(featuresMatrix)) {
        stop("featuresMatrix not of type matrix")
    }
    if (sum(dim(featuresMatrix)) == 2 && is.na(featuresMatrix)) {
        stop("Empty featuresMatrix")
    }

    if (sinuc_or_dinuc == "sinuc") {
        invisible(apply(featuresMatrix, MARGIN = 2, function(x) {
            pwm <- make_sinuc_PWMs(x, add_pseudo_counts = FALSE)
            p1 <- plot_ggheatmap(pwmMat = pwm,
                                    position_labels = position_labels,
                                    savePDFfilename = savePDFfilename)
            base::suppressMessages(print(p1))
        }))
    } else if (sinuc_or_dinuc == "dinuc") {
        invisible(apply(featuresMatrix, MARGIN = 2, function(x) {
            pwm <- make_dinuc_PWMs(x, add_pseudo_counts = FALSE)
            p1 <- plot_ggheatmap(pwmMat = pwm,
                                    position_labels = position_labels,
                                    savePDFfilename = savePDFfilename)
            base::suppressMessages(print(p1))
        }))
    }
}
