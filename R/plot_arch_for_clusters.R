#' @title Plot cluster architectures as sequence logos.
#' 
#' @description Given a collection of FASTA sequences as a DNAStringSet object,
#' and the clusters information, this function plots the architectures for all
#' clusters. If a name for the PDF file is provided, the resulting set of
#' architecture sequence logos are saved as a multi-page PDF.
#'
#' @param seqs Sequences as a \code{\link{Biostrings::DNAStringSet}}.
#' 
#' @param clust_list Clusters as a list of sequence IDs in each cluster.
#' 
#' @param pos_lab Labels for sequence positions, should be of same
#' length as length of the sequences.
#' 
#' @param xt_freq Frequency of x-axis ticks.
#' 
#' @param pdf_width,pdf_height Width and height in inches of the PDF file. 
#' Default values are 11 and 2.
#' 
#' @param pdf_name Specify the PDF filename. 
#' 
#' @param ... Additional arguments passed to the graphic device.
#'
#' @return A list of (ggplot2-based) sequence logo plots is returned. When a 
#' valid file name is specified, the list of plots is also written to the PDF 
#' file (one plot per page).
#'
#' @importFrom grDevices pdf dev.off
#' 
#' @export
plot_arch_for_clusters <- function(seqs,
                                clust_list,
                                pos_lab,
                                xt_freq = 5,
                                pdf_width = 11,
                                pdf_height = 2,
                                pdf_name = "archR_sequence_architectures.pdf",
                                ...){
    ##
    stopifnot(!is.null(clust_list))
    stopifnot(!is.null(seqs))
    if(!is(seqs, "DNAStringSet")){
        stop("Expecting a DNAStringSet object as 'seqs'")
    }
    ##
    nClust <- length(clust_list)
    clust_lens <- unlist(lapply(clust_list, length))
    cumsums_clust_lens <- cumsum(clust_lens)
    clust_names <- sort(as.character(seq_along(clust_list)))
    ##
    clust_starts <- c(1, 1+cumsums_clust_lens[seq_len(nClust-1)])
    clust_ends <- cumsums_clust_lens
    plot_titles <- lapply(seq_along(clust_starts), function(x){
        make_plot_title_str(x, nClust, clust_names[x], clust_lens[x], 
                            clust_starts[x], clust_ends[x])
    })
    ##
    suppressMessages(plot_list <- lapply(seq_along(clust_list), function(x){
        plot_ggseqlogo_of_seqs(seqs=seqs[clust_list[[x]]],
            pos_lab = pos_lab,
            xt_freq = xt_freq,
            title = plot_titles[[x]])
    }))
    
    if(!is.null(pdf_name)){
        grDevices::pdf(file=pdf_name, width=pdf_width, height=pdf_height)
        lapply(plot_list, print)
        dev.off()
    }
    return(plot_list)
}

make_plot_title_str <- function(i, n, name, this_size, st, ed){
    paste0("(", i , "/", n, ") Arch '",
        name, "': ", this_size, " sequences (",  st, "-",  ed, ")")
}

#' @title Plot sequence logo of a collection of sequences
#'
#' @description A wrapper to ggseqlogo plotting. Given a collection of
#' sequences, this function plots the sequence logo.
#'
#' @param seqs Collection of sequences as a Biostrings::DNAStringSet object.
#' 
#' @param pos_lab Numeric vector of labels for sequence positions.
#' This should be the same length as the width of the given sequences.
#' 
#' @param xt_freq Specify the frequency of the x-axis ticks.
#' 
#' @param title The title for the plot.
#' 
#' @param bits_yax Specify 'full' if the information content y-axis limits 
#' should be 0-2 or 'auto' for a suitable limit. The 'auto' setting adjusts 
#' the y-axis limits according to the maximum information content of the 
#' sequence logo. Default is 'auto'.
#' 
#' @importFrom Biostrings width
#' 
#' @export
plot_ggseqlogo_of_seqs <- function(seqs, pos_lab, xt_freq = 5,
                                    title = "Title", bits_yax = "auto"){
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
            method = "bits"
        ) +
        ggplot2::theme_linedraw() +
        ggplot2::theme(axis.text.x = element_text(size = rel(0.9),
                                                    angle = 90,
                                                    hjust = 1),
                        axis.text.y = element_text(size = rel(0.9)),
                        panel.grid = element_blank()
        ) +
        ## Add additional bold tick labels
        ggplot2::scale_x_continuous(breaks = xtick_cal,
                                    labels = pos_lab[xtick_cal],
                                    expand = expansion(mult = c(0, 0))) +
        ggplot2::ggtitle(title)
    if(bits_yax == 'full'){
        foo_p <- foo_p + ggplot2::ylim(0.0, 2.0) 
    }
    message("Plot title:", title)

    foo_p
}
