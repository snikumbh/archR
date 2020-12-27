#' @title Plot cluster architectures as sequence logos.
#' @description Given a collection of FASTA sequences as a DNAStringSet object,
#' and the clusters information, this function plots the architectures for all
#' clusters. If a name for the PDF file is provided, the resulting set of
#' architecture sequence logos are saved as a multi-page PDF.
#'
#' @param tss.seqs_raw Sequences as a DNAStringSet
#' @param list_of_elements Cluster elements as a list
#' @param position_labels Labels for sequence positions, should be of same
#' length as length of the sequences
#' @param xt_freq Frequency of x-axis ticks
#' @param fwidth,fheight Width and height in inches. Default values are 11 and 
#' 2
#' @param PDFfname Specify a filename to be saved as PDF
#'
#' @importFrom grDevices pdf
#' @export
plot_arch_for_clusters <- function(tss.seqs_raw,
                                list_of_elements,
                                position_labels,
                                xt_freq = 5,
                                fwidth = 11,
                                fheight = 2,
                                PDFfname = "archR_sequence_architectures.pdf"){
    if(!is.null(PDFfname)) {
        grDevices::pdf(file=PDFfname, width = fwidth, height = fheight)
    }
    seqs_clusters_as_list <- list_of_elements
    cluster_lengths <- unlist(lapply(seqs_clusters_as_list, length))
    cumsums_of_cluster_lengths <- cumsum(cluster_lengths)
    cluster_names <- sort(as.character(seq_along(seqs_clusters_as_list)))
    for (i in seq_along(seqs_clusters_as_list)) {
        # if(i > 1){
        #     startN <- 1 + cumsums_of_cluster_lengths[i-1]
        # }else{
        #     startN <- 1
        # }
        ##
        startN <- ifelse(i > 1, (1+cumsums_of_cluster_lengths[i-1]), 1)
        endN <- cumsums_of_cluster_lengths[i]
        ##
        plot_title <- paste0("(", i , "/", length(seqs_clusters_as_list), 
                             ") Arch `",
                             cluster_names[i], "': ",
                             length(seqs_clusters_as_list[[i]]),
                             " sequences (",  startN, "-",  endN, ")")
        ##
        foo_p <- plot_ggseqlogo_of_seqs(seqs = 
                                    tss.seqs_raw[ seqs_clusters_as_list[[i]] ],
                                    position_labels = position_labels,
                                    xt_freq = xt_freq,
                                    title = plot_title)
        ##
        suppressMessages(print(foo_p))
    }
    if(!is.null(PDFfname)) {
        dev.off()
    }
}


#' @title Plot sequence logo of a collection of sequences
#'
#' @description A wrapper to ggseqlogo plotting. Given a collection of
#' sequences, this function plots the sequence logo.
#'
#' @param seqs Collection of sequences as a Biostrings::DNAStringSet object.
#' @param position_labels Numeric vector of labels for sequence positions.
#' This should be the same length as the width of the given sequences.
#' @param xt_freq Specify the frequency of the x-axis ticks.
#' @param title The title for the plot.
#' @param bits_yax Specify 'full' if the information content y-axis limits 
#' should be 0-2 or 'auto' for a suitable limit. The 'auto' setting adjusts 
#' the y-axis limits according to the maximum information content of the 
#' sequence logo. Default is 'auto'.
#'
#' @export
plot_ggseqlogo_of_seqs <- function(seqs, position_labels, xt_freq = 5,
                                       title = "Title", bits_yax = "auto"){

    nPos <- length(position_labels)
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
                                    labels = position_labels[xtick_cal],
                                    expand = expansion(mult = c(0, 0))) +
        ggplot2::ggtitle(title)
    if(bits_yax == 'full'){
        foo_p <- foo_p + ggplot2::ylim(0.0, 2.0) 
    }
    message("Plot title:", title)

    foo_p
}
