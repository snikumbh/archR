#' @title Visualize a given (PWM) matrix as a sequence logo.
#'
#' @param pwm_mat Matrix (usually a PWM, but can be any non-normalized matrix) 
#' to be represented as a sequence logo.
#' 
#' @param plot_method For \code{ggseqlogo}; either of 'custom', 'bits', or
#' 'probability'. Default is 'custom'.
#' 
#' @param pos_lab Labels of the positions in the sequences.
#' 
#' @param pdf_fname Name of the file which will be saved as PDF.
#'
#' @return A ggplot2 object so you can simply call \code{print} or \code{save}
#' on it later. If \code{pdf_fname} is given, it is also saved in addition to 
#' returning the ggplot object.
#' 
#' @export
#' 
#' @family visualization functions
#' 
#' @seealso \code{\link{plot_ggheatmap}} for plotting PWMs as heatmaps,
#' \code{\link{plot_ggseqlogo_of_seqs}} for visualizing a collection of 
#' sequences by their sequence logo.
#' 
#' @import ggplot2
#' @import ggseqlogo
plot_ggseqlogo <- function(pwm_mat, plot_method = "custom",
                            pos_lab = NULL, pdf_fname = NULL) {
    ##
    if (!is.matrix(pwmMat)) {
        stop("Expecting a matrix with 4 rows")
    }
    if (sum(dim(pwmMat)) == 2 && is.na(pwmMat)) {
        stop("Empty matrix")
    }
    if (!(nrow(pwmMat) == 4)) {
        stop("Expecting a matrix with 4 rows corresponding to DNA alphabet")
    }
    ##
    if (length(position_labels) < ncol(pwmMat)) {
        stop(paste0("Inadequate position labels supplied",
                    ncol(pwmMat) - length(position_labels)
                    )
            )
    }
    ##
    if (length(position_labels) > ncol(pwmMat)) {
        stop(paste0("Overabundant position labels supplied",
                    length(position_labels) - ncol(pwmMat)
                    )
            )
    }
    ##
    p1 <- ggplot() +
            geom_logo(pwmMat, method = plotMethod, seq_type = "dna") +
            theme_logo() +
            ##
            ggplot2::scale_x_continuous(breaks = seq_len(ncol(pwmMat)),
                                        labels = position_labels,
                                        expand = expansion(mult = c(0, 0))) +
            ##
            ggplot2::theme(axis.text.x = element_text(size = rel(0.5),
                                                angle = 90, hjust = 1),
                            axis.text.y = element_text(size = rel(0.5))) +
            ggplot2::ylim(0.0, 2.0)
    ##
    if (!is.null(savePDFfilename)) {
        if (file.exists(savePDFfilename)) {
            warning("File exists, will overwrite", immediate. = TRUE)
        }
        ggsave(p1, device = "pdf", width = 25, height = 0.5)
    }
    return(p1)
}
