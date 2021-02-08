#' @title Plot cluster architectures as sequence logos.
#' 
#' @description Given a collection of FASTA sequences as a DNAStringSet object,
#' and the clusters information, this function plots the architectures for all
#' clusters. If a name for the PDF file is provided, the resulting set of
#' architecture sequence logos are saved as a multi-page PDF.
#'
#' @param seqs Sequences as a \code{\link[Biostrings]{DNAStringSet}}.
#' 
#' @param clust_list Clusters as a list of sequence IDs in each cluster.
#' 
#' @param pos_lab Labels for sequence positions, should be of same
#' length as that of the sequences. Default value is NULL, when the 
#' positions are labeled from 1 to the length of the sequences.
#' 
#' @param xt_freq Frequency of x-axis ticks.
#' 
#' @param set_titles Specify TRUE if titles are to be written for the plots.
#' With FALSE, there are no titles for the plots. The title for each plot 
#' includes the current cluster number, total number of clusters, start and 
#' end sequence numbers in the collection.
#' 
#' @param pdf_width,pdf_height Width and height in inches of the PDF file. 
#' Default values are 11 and 2.
#' 
#' @param pdf_name Specify the PDF filename. 
#' 
#' @param ... Additional args passed to \code{\link{plot_ggseqlogo_of_seqs}}.
#'
#' @return A list of (ggplot2-based) sequence logo plots is returned. When a 
#' valid file name is specified, the list of plots is also written to the PDF 
#' file (one plot per page).
#'
#' @importFrom grDevices pdf dev.off
#' 
#' @examples 
#' res <- readRDS(system.file("extdata", "example_archRresult.rds", 
#'          package = "archR", mustWork = TRUE))
#' 
#' # Default position labels 1 to length of the sequences.
#' arch_pl <- plot_arch_for_clusters(seqs = seqs_str(res), 
#'                                   clust_list = res$clustSol$clusters,
#'                                   pos_lab = NULL,
#'                                   pdf_name = NULL)
#' arch_pl
#' 
#' # Can also set pos_lab based on biology, e.g., use -50 to 49 denoting 
#' # 50 basepairs upstream and 49 downstream of the transcription start site 
#' # located at position 0.
#' arch_pl <- plot_arch_for_clusters(seqs = seqs_str(res), 
#'                                   clust_list = res$clustSol$clusters,
#'                                   pos_lab = seq(-50,49),
#'                                   pdf_name = NULL)
#' arch_pl
#' 
#' # Plotting architecture sequence logos with probability instead of 
#' # information content
#' arch_pl <- plot_arch_for_clusters(seqs = seqs_str(res), 
#'                                   clust_list = res$clustSol$clusters,
#'                                   pos_lab = seq(-50,49),
#'                                   method = "prob",
#'                                   pdf_name = NULL)
#' arch_pl
#'  
#' @export
plot_arch_for_clusters <- function(seqs,
                                clust_list,
                                pos_lab = NULL,
                                xt_freq = 5,
                                set_titles = TRUE,
                                pdf_width = 11,
                                pdf_height = 2,
                                pdf_name = "archR_sequence_architectures.pdf",
                                ...){
    ##
    stopifnot(!is.null(clust_list))
    stopifnot(!is.null(seqs))
    ## seqs can be a DNAStringSet object or a character vector
    if(is(seqs, "character") && length(seqs) > 0){
        ## OK
    }else if(is(seqs, "DNAStringSet")){
        ## OK
    }else{
        stop("Expecting either a DNAStringSet object or a character vector")
    }
    ##
    if(is.null(pos_lab)){
        pos_lab <- seq_len(Biostrings::width(seqs[1]))
    }
    ##
    plot_titles <- make_plot_titles(clust_list, set_titles)
    suppressMessages(plot_list <- lapply(seq_along(clust_list), function(x){
        plot_ggseqlogo_of_seqs(seqs=seqs[clust_list[[x]]],
            pos_lab = pos_lab,
            xt_freq = xt_freq,
            title = plot_titles[[x]],
            ...)
    }))
    ##
    if(!is.null(pdf_name)){
        grDevices::pdf(file=pdf_name, width=pdf_width, height=pdf_height)
        lapply(plot_list, print)
        dev.off()
    }
    return(plot_list)
}
## =============================================================================

make_plot_titles <- function(clust_list, set_titles){
    if(set_titles){
        nClust <- length(clust_list)
        clust_lens <- unlist(lapply(clust_list, length))
        cumsums_clust_lens <- cumsum(clust_lens)
        clust_names <- sort(as.character(seq_along(clust_list)))
        ##
        clust_starts <- c(1, 1+cumsums_clust_lens[seq_len(nClust-1)])
        clust_ends <- cumsums_clust_lens
        ##
        pl_titles <- lapply(seq_along(clust_starts), function(x){
            make_plot_title_str(x, nClust, clust_names[x], clust_lens[x], 
                clust_starts[x], clust_ends[x])
        })
    }else{
        pl_titles <- lapply(seq_along(clust_starts), function(x) NULL)
    }
    pl_titles
}
## =============================================================================

make_plot_title_str <- function(i, n, name, this_size, st, ed){
    paste0("(", i , "/", n, ") Arch '",
        name, "': ", this_size, " sequences (",  st, "-",  ed, ")")
}
## =============================================================================

#' @title Plot sequence logo of a collection of sequences
#'
#' @description A wrapper to ggseqlogo plotting. Given a collection of
#' sequences, this function plots the sequence logo.
#'
#' @param seqs Collection of sequences as a 
#' \code{\link[Biostrings]{DNAStringSet}} object.
#' 
#' @param pos_lab Labels for sequence positions, should be of same
#' length as that of the sequences. Default value is NULL, when the 
#' positions are labeled from 1 to the length of the sequences.
#' 
#' @param xt_freq Specify the frequency of the x-axis ticks.
#' 
#' @param method Specify either 'bits' for information content or 
#' 'prob' for probability.
#' 
#' @param title The title for the plot.
#' 
#' @param bits_yax Specify 'full' if the information content y-axis limits 
#' should be 0-2 or 'auto' for a suitable limit. The 'auto' setting adjusts 
#' the y-axis limits according to the maximum information content of the 
#' sequence logo. Default is 'auto'.
#' 
#' @return A sequence logo plot of the given DNA sequences.
#' 
#' @seealso \code{\link{plot_arch_for_clusters}} for obtaining multiple 
#' sequence logo plots as a list.
#' 
#' @importFrom Biostrings width
#' 
#' @examples 
#' res <- readRDS(system.file("extdata", "example_archRresult.rds", 
#'          package = "archR", mustWork = TRUE))
#' 
#' # Default, using information content on y-axis
#' pl <- plot_ggseqlogo_of_seqs(seqs = seqs_str(res, iter=1, cl=1),
#'                              pos_lab = seq_len(1:100))
#'                              
#' # Using probability instead of information content
#' pl <- plot_ggseqlogo_of_seqs(seqs = seqs_str(res, iter=1, cl=1),
#'                              pos_lab = seq_len(1:100), 
#'                              method = "prob")
#'                              
#' @export
plot_ggseqlogo_of_seqs <- function(seqs, pos_lab = NULL, xt_freq = 5, 
                                    method = "bits", title = "Title", 
                                    bits_yax = "auto"){
    ##
    if(is.null(pos_lab)){
        pos_lab <- seq_len(Biostrings::width(seqs[1]))
    }
    ##
    if(xt_freq > Biostrings::width(seqs[1])){
        xt_freq <- 5
    }
    ##
    nPos <- length(pos_lab)
    xtick_cal <- seq(0, nPos, by = xt_freq)
    xtick_cal[1] <- 1
    xtick_cal[length(xtick_cal)] <- nPos

    foo_p <-
        ggseqlogo::ggseqlogo(
            as.character(seqs),
            seq_type = "dna",
            method = method
        ) +
        ggplot2::theme_linedraw() +
        ggplot2::theme(axis.text.x = element_text(size = rel(0.9),
                                        angle = 90, hjust = 1, vjust=0.5),
                        axis.text.y = element_text(size = rel(0.9)),
                        panel.grid = element_blank()
        ) +
        ## Add additional bold tick labels
        ggplot2::scale_x_continuous(breaks = xtick_cal,
                                    labels = pos_lab[xtick_cal],
                                    expand = expansion(mult = c(0, 0)))
    ##
    if(!is.null(title)) foo_p <- foo_p + ggplot2::ggtitle(title)
    if(bits_yax == 'full') foo_p <- foo_p + ggplot2::ylim(0.0, 2.0) 
    ##
    .msg_pstr("Plot title:", title, flg=TRUE)
    foo_p
}
## =============================================================================
