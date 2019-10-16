#' @title Getter function for labels of x- or y-ticks.
#'
#' @description A convenience function to provide the right set of x-tick labels
#'  for the image matrix representing the sequences.
#'
#'
#' @return A label vector that can be directly used for setting the tickmarks.
#'
get_tick_labels <- function(givenVec, use_asis = TRUE, at_every) {
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


seqs_to_df <- function(seqs, position_labels, annClusters = NULL,
                       sinuc_or_dinuc = "sinuc"){
    # initialize data frames beforehand?
    # temp_df <- data.frame(nucleotides = character(length(position_labels)),
    #                       seq_id = numeric(length(position_labels)),
    #                       positions = numeric(length(position_labels)))

    if(sinuc_or_dinuc == "sinuc"){
        ## handling single nucleotides
        plot_df_list <- lapply(1:length(seqs),
                          function(x) {
                            # given one hot encoded seq to str on nucleotides
                            str_seq <- seqs[x]
                            nucleotides <- unlist(strsplit(str_seq,
                                                           split = NULL))
                            temp_df <- as.data.frame(nucleotides,
                                                     stringsAsFactors = FALSE)
                            ## This is before adding other columns, because
                            ## The last state of the variable is returned
                            if(!is.null(annClusters)){
                              temp_df <- tibble::add_column(temp_df,
                                                            annClusters =
                                                              annClusters[x])
                            }
                            ##
                            temp_df <- tibble::add_column(temp_df,
                                                          seq_id = x,
                                                          positions =
                                                            position_labels)
                          }
                          )
    } else if (sinuc_or_dinuc == "dinuc"){
        ## handling dinucleotides
        plot_df_list <-
          lapply(seq_along(seqs),
                 function(x){
                    x_split <- unlist(strsplit(seqs[x],
                                               split = NULL))
                    this_seq_len <- length(x_split)
                    nucleotides <- unlist(lapply(seq_len(this_seq_len - 1),
                                                 function (y){
                                                   paste0(x_split[y],
                                                          x_split[y+1])
                                                 }
                                                )
                                          )
                    temp_df <- as.data.frame(nucleotides,
                                             stringsAsFactors = FALSE)
                    if(!is.null(annClusters)) {
                       temp_df <- tibble::add_column(temp_df, annClusters =
                                                    annClusters[x])
                    }
                    ##
                    position_labels_dinuc <-
                      position_labels[seq_len(length(position_labels) - 1)]
                    temp_df <- tibble::add_column(temp_df,
                                         seq_id = x,
                                         positions = position_labels_dinuc)
                  }
                )
    }
    plot_df <- data.table::rbindlist(plot_df_list)
    return(plot_df)
}


get_seq_cluster_breaks <- function(annClusters){

  clusterLevels <- levels(as.factor(annClusters))
  interimBreaks <- vapply(seq_along(clusterLevels), function(x){
                          max(which(annClusters == clusterLevels[x]))
                          }, numeric(1))
  seqClustBreaks <- c(1, interimBreaks)

  return(seqClustBreaks)
}

#' @title samarth
#' @description samarth
#' @param k number of colors
#' @importFrom grDevices col2rgb rgb
.distinctColorPalette <- function(k) {
    # set.seed(123)
    ColorSpace <- t(unique(col2rgb(scales::hue_pal(l = 60:100)(2e3))))
    km <- kmeans(ColorSpace, k, iter.max = 20)
    colors <- rgb(round(km$centers), maxColorValue = 255)
    return(colors)
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
#' percentage N-th position. For every N-th position, simply specify the
#' absolute positive value. For every percentage N-th position, specify 0.01 for
#'  10\%.
#' @param plot.title The title of the plot.
#' @param saveFilename Name of the file which will be saved as PDF.
#' @param fileType give file type.
#' @param sinuc_or_dinuc Specify single- or dinucleotide profiles.
#' @param annClusters Give annotation clusters.
#' @param choose_colors OK.
#' @param verbose Default \code{0} which will not print any messages, or can be
#' set to \code{1} which will print messages.
#'
#' @return nothing.
#'
#' @export
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom rlang .data
represent_matrix_of_acgt_with_ggplot2 <- function(givenMat,
                                     position_labels,
                                     sinuc_or_dinuc = "sinuc",
                                     choose_colors = c("A" = "darkgreen",
                                                       "C" = "blue",
                                                       "G" = "orange",
                                                       "T" = "red"),
                                     annClusters = NULL,
                                     mark_seq_at_every = as.integer(
                                       ncol(givenMat) * 0.1),
                                     # percentage value, 10 as default
                                     mark_pos_at_every = 10,
                                     plot.title = "DNA Sequences",
                                     saveFilename = "DNA_sequences_as_matrix",
                                     fileType = "png",
                                     verbose = 0) {
  #
  # Returns: No return value
  # nPositions <- nrow(givenMat)/4
  nSeqs <- ncol(givenMat)
  # mat is now of dim 4 x positions
  if (verbose > 0) {
    cat("Plotting\n")
  }
  # iupac_colors <- c("A" = "darkgreen",
  #                   "C" = "blue",
  #                   "G" = "orange",
  #                   "T" = "red")

  plot_df <- seqs_to_df(givenMat,
                        position_labels = position_labels,
                        annClusters = annClusters,
                        sinuc_or_dinuc = sinuc_or_dinuc)

  #
  # Position Breaks on axis
  position_levels <- as.numeric(levels(as.factor(plot_df$positions)))
  position_breaks <- position_levels[which(position_levels%%5 == 0)]
  if (! (min(position_levels) %in% position_breaks)){
    position_breaks <- c(1, position_breaks)

  }

  # Seqeunce Breaks on axis
  sequence_levels <- as.numeric(levels(as.factor(plot_df$seq_id)))
  sequence_breaks <- which(sequence_levels%%100 == 0)
  if (! (1 %in% sequence_breaks)){
    sequence_breaks <- c(1, sequence_breaks)

  }
  ###

  pMatrix <- ggplot2::ggplot(data = plot_df, mapping = aes(
    x = .data$positions,
    # Here, 'positions' is the column_name, see previous statement.
    # Do not change it to position_labels
    y = .data$seq_id,
    fill = .data$nucleotides
  )) +
    ggplot2::geom_tile() +
    ggplot2::theme_bw() +
    ggplot2::xlab(label = "Positions") +
    ggplot2::ylab(label = 'Sequences') +
    ggplot2::scale_y_continuous(expand = c(0.01, 0.01),
                                trans = "reverse",
                                position = "left",
                                breaks = sequence_breaks) +
    #
    ggplot2::scale_x_continuous(expand = c(0.01,0.01),
                                position = "bottom",
                                breaks = position_breaks) +
    #
    ggplot2::scale_fill_manual(values = choose_colors) +
    #ggplot2::coord_fixed(ratio = 0.05) +
    ggplot2::labs(title = plot.title) +
    ggplot2::theme(
      legend.position = "left",
      legend.title = element_blank(),
      axis.text.x = element_text(size=rel(0.9), angle = 90, hjust = 1),
      axis.text.y = element_text(size=rel(0.9))
    )
  ###
  if(!is.null(annClusters)){
      pCluster <- ggplot2::ggplot(data = plot_df,
                                  mapping = aes(x = 0,
                                  y = .data$seq_id,
                                  fill = .data$annClusters)) +
        ggplot2::geom_tile() +
        ggplot2::theme_bw() +
        ggplot2::ylab(label = 'Clusters') +
        ggplot2::xlab(label = NULL) +
        ggplot2::theme(legend.position = 'none') +
        ggplot2::scale_y_continuous(expand = c(0.01, 0.01),
                                    trans = "reverse",
                                    position = "right",
                                    breaks =
                                      get_seq_cluster_breaks(annClusters)) +
        ggplot2::scale_x_continuous(expand = c(0.0,0.0),
                                    position = "bottom",
                                    labels = NULL, breaks = NULL) +
        ggplot2::scale_fill_manual(values = .distinctColorPalette(length(
                                        levels(as.factor(annClusters))
                                  )))
        # ggplot2::scale_fill_manual(values = sample(get_palette(
        #                                         palette = "Paired",
        #                                        length(levels(
        #                                          as.factor(annClusters))))))

      ##
      final_p <- ggpubr::ggarrange(pMatrix, pCluster, ncol = 2, heights = 1:1,
                                   widths = c(8,1.5), align = 'h')
  } else{
      final_p <- pMatrix
  }
  ###
  if (!is.null(saveFilename)) {
    if(fileType == "pdf"){
      ggplot2::ggsave(paste0(saveFilename, ".pdf"),
                      final_p, device = "pdf", width = 6.0, height = 11.69)
    } else if (fileType == "png"){
      ggplot2::ggsave(paste0(saveFilename, ".png"),
                      final_p, device = "png", width = 6.0, height = 11.69)
    }
  }else{
      print(final_p)

  }

}
