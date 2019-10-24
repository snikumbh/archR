#' @title Plot a given matrix as a heatmap using ggplot2.
#'
#' @description The given matrix is plotted as a heatmap using \code{ggplot2}'s
#' \code{geom_tile}.
#'
#' @param pwmMat Matrix (usually a PWM, but can be non-normalized/any matrix)
#' to be represented as a heatmap.
#' @param position_labels Labels of the positions in the sequences.
#' @param savePDFfilename Name of the file which will be saved as PDF.
#'
#' @return A ggplot object so you can simply call \code{print} or \code{save}
#' on it later. If \code{savePDFfilename} is given, it is also saved and the
#' \code{ggplot} object returned.
#' @export
#' @family visualization functions
#' @seealso \code{\link{plot_ggseqlogo}} for plotting PWMs as sequence logos
#' @importFrom reshape2 melt
#' @importFrom tibble add_column
#' @import ggplot2
#' @import ggseqlogo
#'
plot_ggheatmap <- function(pwmMat, position_labels = NULL,
                            savePDFfilename = NULL) {
    if (!is.matrix(pwmMat)) {
        stop("Expecting a matrix with 4 rows")
    }
    if (sum(dim(pwmMat)) == 2 && is.na(pwmMat)) {
        stop("Empty matrix")
    }
    if (!(nrow(pwmMat) == 4)) {
        stop("Expecting a matrix with 4 rows corresponding to DNA chars ",
                "'A', 'C', 'G', 'T'")
    }
    ##
    if (length(position_labels) < ncol(pwmMat)) {
        stop(paste0(
            "Inadequate position labels supplied",
            ncol(pwmMat) - length(position_labels)
        ))
    }
    ##
    if (length(position_labels) > ncol(pwmMat)) {
        stop(paste0(
            "Overabundant position labels supplied",
            length(position_labels) - ncol(pwmMat)
        ))
    }
    ##
    ## Convert pwmMat to df, heatmap by ggplot-way
    pwmMat_df <- as.data.frame(pwmMat)

    pwmMat_df <- add_column(pwmMat_df, Nucleotides = rownames(pwmMat_df))

    colnames(pwmMat_df) <- c(position_labels, "Nucleotides")
    ##
    pwmMat_df_for_ggheatmap <- melt(pwmMat_df, id.vars = c("Nucleotides"),
                                    variable.name = "positions")
    ##
    p1 <- ggplot2::ggplot(data = pwmMat_df_for_ggheatmap, mapping = aes(
        x = positions,
        ## Here, 'positions' is the column_name, see previous statement.
        ## Do not change it to position_labels
        y = Nucleotides,
        fill = value
    )) +
        ggplot2::geom_tile() +
        ggplot2::theme_bw() +
        ggplot2::xlab(label = "Positions") +
        ggplot2::scale_fill_gradient2(
            name = "",
            low = "white",
            mid = "white",
            high = "#012345"
        ) +
        ggplot2::coord_fixed(ratio = 2.0, clip = "on") +
        ggplot2::theme(
            legend.position = "top",
            legend.justification = "center",
            axis.text.x = element_text(size = rel(0.5), angle = 90, hjust = 1)
        )

    if (!is.null(savePDFfilename)) {
        if (file.exists(savePDFfilename)) {
            warning("File exists, will overwrite")
        }
        ggplot2::ggsave(p1, device = "pdf", width = 20, height = 2.5)
    }
    return(p1)
}

