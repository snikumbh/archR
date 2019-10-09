#' @title Getter function for labels of x- or y-ticks.
#'
#' @description A convenience function to provide the right set of x-tick labels
#'  for the image matrix representing the sequences.
#'
#' @param givenVec The vector of values for the x- or y-ticks.
#' @param use_asis Use the vector \emph{as is}, i.e., all given values will
#' be shown as ticks.
#' @param at_every Show tickmarks only at positions specified.
#'
#' @return A label vector that can be directly used for setting the tickmarks.
#'
get_tick_labels <- function(givenVec, use_asis = T, at_every) {
  givenLength <- length(givenVec)
  label_vec <- rep("", givenLength)
  if (use_asis) {
    indices <- seq(0, givenLength / at_every) * at_every + 1
    indices[length(indices)] <- length(givenVec)
    values <- givenVec[indices]
    label_vec[indices] <- values
  } else {
    indices <- seq(givenLength / at_every) * at_every
    values <- seq(givenLength / at_every) * at_every
    label_vec[indices] <- values
  }
  return_val <- list(label_ind = 1:givenLength, label_val = label_vec)
  return(return_val)
}


#' @title Represent a Matrix of Sequences as an Image.
#'
#' @description The data matrix is expected to have features alogn the rows and
#' sequences along the columns.
#'
#' @param givenMat The data matrix describing the set of sequences to be
#' visualized as an image. The current version is applicable only for DNA
#' sequences.
#' @param position_labels Labels of the positions in the sequences.
#' @param mark_seq_at_every Tickmarks for sequences, i.e., these are y-ticks.
#' @param mark_pos_at_every Tickmarks for sequence positions, i.e., these are
#' x-ticks. The tickmarks can be specified for every N-th position or every
#' percentage N-th position. For every N-th position, simply specify the a
#' bsolute positive value. For every percentage N-th position, specify 0.01 for
#' 10\%.
#' @param plot.title The title of the plot.
#' @param savePDFfilename Name of the file which will be saved as PDF.
#' @param verbose Default \code{0} which will not print any messages, or can be
#' set to \code{1} which will print messages.
#'
#' @return nothing.
#'
#' @export
#' @importFrom grDevices dev.off pdf colorRampPalette
#'
represent_matrix_of_acgt <- function(givenMat,
                                     position_labels,
                                     mark_seq_at_every =
                                       as.integer(ncol(givenMat) * 0.1),
                                     # percentage value, 10 as default
                                     mark_pos_at_every = 10,
                                     plot.title = "DNA Sequences",
                                     savePDFfilename = NULL,
                                     verbose = 0) {
  # givenMat is expected as a matrix of #Features x #Sequences
  # WARNING: This seems to not work when pseudo-counts are added to the data
  # matrix from generate_seq_info_matrix function
  #
  # Returns: No return value
  # nPositions <- nrow(givenMat)/4
  nSeqs <- ncol(givenMat)
  sinuc <- c("A", "C", "G", "T")
  iupac_colors <- c("darkgreen", "blue", "orange", "red")
  # mat is now of dim 4 x positions
  plot_mat <- apply(givenMat, MARGIN = 2, function(x) {
    # print(dim(as.matrix(x)))
    this_givenMatrix <- matrix(x, nrow = length(sinuc), byrow = TRUE)
    rownames(this_givenMatrix) <- sinuc
    new_mat <- matrix(rep(0, ncol(this_givenMatrix)),
                      ncol = ncol(this_givenMatrix))
    for (r in 1:4) {
      temp_mat <- this_givenMatrix
      temp_mat[-r, ] <- 0
      new_mat[, which(temp_mat[r, ] == 1) ] <- r # sinuc[r]
    }
    return(new_mat)
  })
  # To get back the sequence
  # paste0(rs,collapse='')
  #
  if (verbose > 0) {
    cat("Creating a color palette based IUPAC colors for nucleotides\n")
  }
  acgt_palette <- grDevices::colorRampPalette(iupac_colors)(n = 4)
  # given matrix is #Features x #Sequences
  # we plot #Sequences x #Features
  plot_mat <- t(plot_mat)
  #
  if (!is.null(savePDFfilename)) {
    pdf(savePDFfilename)
  }
  # heatmap.2 solution is slow (available in `
  # represent_matrix_of_acgt_heatmap2.R`),
  # we instead use heatmap3 which is fast
  #
  if (verbose > 0) {
    cat("Plotting\n")
  }
  # heatmap3 solution
  seq_labels <- get_tick_labels(1:nSeqs, use_asis = F,
                                at_every = mark_seq_at_every)
  pos_labels <- get_tick_labels(position_labels, use_asis = T,
                                at_every = mark_pos_at_every)
  heatmap3::heatmap3(plot_mat,
    revC = T, # set this to T to plot starting from the top
    Colv = NA,
    Rowv = NA,
    col = acgt_palette,
    main = plot.title,
    margins = c(2.5, 3),
    scale = "none",
    # positions along x-axis
    xlab = paste0("positions"),
    labCol = pos_labels$label_val,
    # sequences along y-axis
    ylab = paste0("sequences"),
    labRow = seq_labels$label_val,
    # xtick = T,
    # add.expr = function(pos_labels=pos_labels, seq_labels=seq_labels) {
    #                        axis(side=1, at = pos_labels$label_val,
    #                        labels = pos_labels$label_val, tick = T)
    #                        #axis(side=2, at = , labels = seq_labels, tick = T)
    #                       },
    legendfun = function() heatmap3::showLegend(
        legend = sinuc,
        lwd = NA,
        horiz = F,
        seg.len = 0.1,
        fill = iupac_colors,
        title = "Nucleotides" # Legend/color key title
      ),
    useRaster = TRUE
  )
  # axis(side=1, at = pos_labels$label_val,
  # labels = pos_labels$label_val, tick = T)
  # axis(side=2, at = seq_labels$label_ind,
  # labels = seq_labels$label_val, tick = T)
  if (!is.null(savePDFfilename)) {
    grDevices::dev.off()
  }
}
