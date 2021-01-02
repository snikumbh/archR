# @title
# Convert DNA sequences from DNAStringSet object to a numeric matrix
#
# @description This function converts a DNAStringSet object into a numeric
# matrix
# @param seqs The sequences to converted to matrix
# @param position_labels Labels for the axis denoting sequence positions
# @param annClusters For later, when clusters are to be annotated
# @param sinuc_or_dinuc For later, maybe
.seqs_to_mat <- function(seqs, position_labels, annClusters = NULL,
                            sinuc_or_dinuc = "sinuc") {
    nSeqs <- length(seqs)
    nPos <- length(position_labels)
    ## handling single- and di-nucleotides separately
    if (sinuc_or_dinuc == "sinuc") {
        nuc_list <- unlist(lapply(seq_along(seqs),
                                function(x) {
                                    str_seq <- seqs[x]
                                    nucleotides <-
                                        unlist(strsplit(str_seq,
                                                split = NULL))
                                }))
        nuc_list[which(nuc_list == "A")] <- 1
        nuc_list[which(nuc_list == "C")] <- 2
        nuc_list[which(nuc_list == "G")] <- 3
        nuc_list[which(nuc_list == "T")] <- 4
        nuc_list_num <- as.numeric(nuc_list)
        nuc_mat <- matrix(nuc_list_num, byrow = TRUE, ncol = nPos, nrow = nSeqs)
    } else if (sinuc_or_dinuc == "dinuc") {
        ## handling dinucleotides
        ## not needed any more?!
        stop("Nothing for dinuc here. Unrelated!")
    }
    return(nuc_mat)
}
## =============================================================================


#' @title
#' Visualize raw DNA sequences as an image
#'
#' @description This function plots the collection of sequences as an image
#' matrix.
#'
#' @param rawSeqs The sequences as a DNAStringSet object.
#' @param position_labels The labels to be used for the sequence positions.
#' Default: Sequence positions are labeled from 1 to the length of the
#' sequences.
#' @param xt_freq The x-axis tick frequency.
#' @param yt_freq The y-axis tick frequency.
#' @param col A vector of four colors used for the DNA bases A, C, G, and T (in
#' that order)
#' @param savefilename Specify the filename (with extension) for saving the
#' plot to disk.
#' @param filetype Specify the file type, namely PNG, JPEG, TIFF.
#' @param fwidth Specify the width for the plot. This depends on the length of
#' sequences.
#' @param fheight Specify the height for the plot. This depends on the number of
#' sequences.
#' @param funits Specify the units in which the height and width are given.
#'
#' @importFrom Biostrings width
#'
#' @return Nothing returned to the R interpreter.
#' @family visualization functions
#' @importFrom grDevices png dev.off
#' @importFrom graphics axis image
#' @export
viz_seqs_as_acgt_mat_from_seqs <- function(rawSeqs, position_labels = NULL,
                                    xt_freq = 5, yt_freq = 100,
                                    col = c("darkgreen", "blue",
                                            "orange", "red"),
                                    savefilename = NULL,
                                    filetype = "PNG",
                                    fwidth = 450, fheight = 900, funits = "px"
                                    ) {
    if(is.null(position_labels)){
        position_labels <- seq_len(Biostrings::width(rawSeqs[1]))
    }

    nSeqs <- length(rawSeqs)
    nPos <- length(position_labels)
    seq_mat <- .seqs_to_mat(seqs = rawSeqs, position_labels = position_labels)
    seq_mat <- seq_mat[rev(seq_len(nSeqs)),]

    if(!is.null(savefilename)){
        if(filetype == "PNG" || filetype == "png"){
            grDevices::png(filename = savefilename,
                width = fwidth, height = fheight, units = funits,
                bg = "white")
        }else{
            ## TODO
        }
    }

    xtick_cal <- seq(0, nPos, by = xt_freq)
    xtick_cal[1] <- 1
    ##
    ytick_names <- rev(seq(yt_freq, nSeqs, by = yt_freq))
    ytick_loc <- 1 + nSeqs - c(rev(seq(yt_freq, nSeqs, by = yt_freq)))

    graphics::image(x = seq_len(nPos), y = seq_len(nSeqs),
                z = t(seq_mat),
                col = col,
                useRaster = TRUE,
                ylab = paste0("Sequences (n = ", nSeqs, ")"),
                xlab = "Positions",
                axes = FALSE)
    axis(side = 1, at = xtick_cal, labels = position_labels[xtick_cal], las = 2)
    axis(side = 2, at = ytick_loc, labels = ytick_names, las = 2)

    if(!is.null(savefilename)){
        dev.off()
    }

}
