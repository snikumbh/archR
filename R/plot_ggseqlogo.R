#' @title Visualize a given (PWM) matrix as a sequence logo.
#'
#' @param pwm_mat Matrix (usually a PWM, but can be any non-normalized matrix) 
#' to be represented as a sequence logo.
#' 
#' @param method For \code{ggseqlogo}; either of 'custom', 'bits', or
#' 'probability'. Default is 'custom'.
#' 
#' @param pos_lab Labels of the positions in the sequences.
#' 
#' @param pdf_name Name of the file which will be saved as PDF.
#'
#' @return A ggplot2 object so you can simply call \code{print} or \code{save}
#' on it later. If \code{pdf_name} is given, it is also saved in addition to 
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
plot_ggseqlogo <- function(pwm_mat, method = "custom", pos_lab = NULL, 
                            pdf_name = NULL) {
    ##
    if (!is.matrix(pwm_mat)) {
        stop("Expecting a matrix with 4 rows")
    }
    if (sum(dim(pwm_mat)) == 2 && is.na(pwm_mat)) {
        stop("Empty matrix")
    }
    if (!(nrow(pwm_mat) == 4)) {
        stop("Expecting a matrix with 4 rows corresponding to DNA alphabet")
    }
    ##
    if (length(pos_lab) < ncol(pwm_mat)) {
        stop(paste0("Inadequate position labels supplied",
                    ncol(pwm_mat) - length(pos_lab)
                    ))
    }
    ##
    if (length(pos_lab) > ncol(pwm_mat)) {
        stop(paste0("Overabundant position labels supplied",
                    length(pos_lab) - ncol(pwm_mat)
                    )
            )
    }
    ##
    p1 <- ggplot() +
            geom_logo(pwm_mat, method = method, seq_type = "dna") +
            theme_logo() +
            ##
            ggplot2::scale_x_continuous(breaks = seq_len(ncol(pwm_mat)),
                                        labels = pos_lab,
                                        expand = expansion(mult = c(0, 0))) +
            ##
            ggplot2::theme(axis.text.x = element_text(size = rel(0.5),
                                                angle = 90, hjust = 1),
                            axis.text.y = element_text(size = rel(0.5)))
    ##
    if (!is.null(pdf_name)) {
        if (file.exists(pdf_name)) {
            warning("File exists, will overwrite", immediate. = TRUE)
        }
        ggsave(p1, device = "pdf", width = 25, height = 0.5)
    }
    return(p1)
}

check_vars <- function(pwm_mat, pos_lab){
    if (!is.matrix(pwm_mat)) {
        stop("Expecting a matrix with 4 rows")
    }
    if (sum(dim(pwm_mat)) == 2 && is.na(pwm_mat)) {
        stop("Empty matrix")
    }
    if (!(nrow(pwm_mat) == 4)) {
        stop("Expecting a matrix with 4 rows corresponding to DNA chars ",
            "'A', 'C', 'G', 'T'")
    }
    ##
    if (length(pos_lab) < ncol(pwm_mat)) {
        stop(paste0("Inadequate position labels supplied", 
            ncol(pwm_mat) - length(pos_lab)
        ))
    }
    ##
    if (length(pos_lab) > ncol(pwm_mat)) {
        stop(paste0("Overabundant position labels supplied",
            length(pos_lab) - ncol(pwm_mat)
        ))
    }
}
