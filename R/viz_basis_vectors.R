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
#' @param method For \code{ggseqlogo} -- either of "custom", "bits", or
#' "probability". Default is "custom".
#' @param pos_lab Labels for sequence positions, should be of same
#' length as that of the sequences. Default value is NULL, when the 
#' positions are labeled from 1 to the length of the sequences.
#' @param add_pseudo_counts Boolean, taking values TRUE/T or FALSE/F, default
#' set to FALSE. Setting it to TRUE will enable adding pseudo-counts to the
#' features matrix.
#' @param pdf_name Name of the file which will be saved as PDF
#' (also provide the extension).
#' @param sinuc_or_dinuc "sinuc" or "dinuc"
#'
#' @return nothing
#'
#' @export
#' @family visualization functions
#'
#' @import ggplot2
#' @import ggseqlogo
#' 
#' @examples 
#' res <- readRDS(system.file("extdata", "example_archRresult.rds", 
#'          package = "archR", mustWork = TRUE))
#' 
#' viz_basis_vectors_in_combined_heatmaps_seqlogos(
#'                   feat_mat = get_clBasVec_m(res,iter=1))
#' 
#' 
viz_basis_vectors_in_combined_heatmaps_seqlogos <-
    function(feat_mat, method = "custom", pos_lab = NULL, 
                add_pseudo_counts = FALSE, pdf_name = NULL,
                sinuc_or_dinuc = "sinuc") {
    check_ggpubr()
    check_vars2(feat_mat)
    ##
    if(is.null(pos_lab)){
        pos_lab <- set_default_pos_lab(feat_mat, sinuc_or_dinuc)
    }
    ##
    invisible(apply(feat_mat, MARGIN = 2, function(x) {
    if (sinuc_or_dinuc == "dinuc") {
        pwm <- make_dinuc_PWMs(x, add_pseudo_counts = FALSE)
    } else if (sinuc_or_dinuc == "sinuc") {
        pwm <- make_sinuc_PWMs(x, add_pseudo_counts = FALSE)
    }
    ## Heatmap on top
    p1 <- plot_ggheatmap(pwm_mat = pwm, pos_lab = pos_lab,
                    pdf_name = pdf_name)
    p1 <- p1 + theme(plot.margin = grid::unit(c(0, 0, 0, 0), "mm"))
    ## Seqlogo below
    p2 <- plot_ggseqlogo(pwm_mat = pwm, method = method,
        pos_lab = pos_lab, pdf_name = pdf_name)
    ## Make adjustments for alignment
    p2 <- p2 + theme(plot.margin = grid::unit(c(0, 0, 0, 0), "mm"))
    final_p <- ggpubr::ggarrange(p1, p2, ncol = 1, heights = c(1.2,1),
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



#' @describeIn viz_basis_vectors_in_combined_heatmaps_seqlogos Visualize the 
#' NMF basis vectors as a sequence logo
#' 
#' @return nothing.
#'
#' @examples 
#' res <- readRDS(system.file("extdata", "example_archRresult.rds", 
#'          package = "archR", mustWork = TRUE))
#' 
#' viz_basis_vectors_as_seqlogo(feat_mat = get_clBasVec_m(res,iter=1))
#'          
#'
#' @export
viz_basis_vectors_as_seqlogo <- function(feat_mat,
                                            method = "custom",
                                            pos_lab = NULL,
                                            add_pseudo_counts = FALSE,
                                            pdf_name = NULL,
                                            sinuc_or_dinuc = "sinuc") {
    ## Visualize all basis factors (expected as columns of the given features
    ## matrix) as seqlogos
    check_vars2(feat_mat)
    ##
    if(is.null(pos_lab)){
        pos_lab <- set_default_pos_lab(feat_mat, sinuc_or_dinuc)
    }
    ##
    invisible(apply(feat_mat, MARGIN = 2, function(x) {
        if (sinuc_or_dinuc == "dinuc") {
            pwm <- make_dinuc_PWMs(x, add_pseudo_counts = FALSE)
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

#' @describeIn viz_basis_vectors_in_combined_heatmaps_seqlogos Visualize the 
#' NMF basis vectors as a heatmap
#' 
#' @return nothing.
#'
#' @examples 
#' res <- readRDS(system.file("extdata", "example_archRresult.rds", 
#'          package = "archR", mustWork = TRUE))
#' 
#' viz_basis_vectors_as_seqlogo(feat_mat = get_clBasVec_m(res,iter=1))
#'          
#'
#' @export
viz_basis_vectors_as_heatmap <- function(feat_mat, pos_lab = NULL,
                                        add_pseudo_counts = FALSE,
                                        pdf_name = NULL,
                                        sinuc_or_dinuc = "sinuc") {
    # Visualize all basis factors (expected as columns of the given features
    # matrix) as heatmaps
    ##
    check_vars2(feat_mat)
    ##
    if(is.null(pos_lab)){
        pos_lab <- set_default_pos_lab(feat_mat, sinuc_or_dinuc)
    }
    ##
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
## =============================================================================

set_default_pos_lab <- function(feat_mat, sinuc_or_dinuc){
    pos_lab <- NULL
    if(sinuc_or_dinuc == "sinuc"){
        pos_lab <- seq_len(nrow(feat_mat)/4)
    }
    if(sinuc_or_dinuc == "dinuc"){
        pos_lab <- seq_len(nrow(feat_mat)/16)
    }
    pos_lab
}

check_ggpubr <- function(){
    if(!requireNamespace("ggpubr", quietly = TRUE)){
        stop("Please install the R package 'ggpubr' to use this function")
    }
}
## =============================================================================

check_vars2 <- function(feat_mat){
    if (!is.matrix(feat_mat)) {
        stop("feat_mat not of type matrix")
    }
    if (sum(dim(feat_mat)) == 2 && is.na(feat_mat)) {
        stop("Empty feat_mat")
    }
}
## =============================================================================
