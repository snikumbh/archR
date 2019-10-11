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
#'
#' @importFrom dplyr mutate
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
    # if (!(nrow(pwmMat) == 16)) {
        stop("Expecting a matrix with 4 rows corresponding to DNA alphabet")
        # stop("Expecting a matrix with 4 or 8 rows corresponding to
        #      DNA alphabet or dinucleotide prfiles respectively")
    # }
  }
  #
  if (length(position_labels) < ncol(pwmMat)) {
    stop(paste0(
      "Inadequate position labels supplied",
      ncol(pwmMat) - length(position_labels)
    ))
  }
  #
  if (length(position_labels) > ncol(pwmMat)) {
    stop(paste0(
      "Overabundant position labels supplied",
      length(position_labels) - ncol(pwmMat)
    ))
  }
  #
  # Convert pwmMat to df, heatmap by ggplot-way
  pwmMat_df <- as.data.frame(pwmMat)

  pwmMat_df <- add_column(pwmMat_df, Nucleotides = rownames(pwmMat_df))

  # pwmMat_df <- mutate(pwmMat_df, Nucleotides = rownames(pwmMat_df))
  # colnames(pwmMat_df) <- c(position_labels, "Nucleotides")
  # WHY THE HELL IS THIS THIS WAY?
  colnames(pwmMat_df) <- c(position_labels, "Nucleotides")
  #
  pwmMat_df_for_ggheatmap <- melt(pwmMat_df, id.vars = c("Nucleotides"),
                                  variable.name = "positions")
  #
  p1 <- ggplot2::ggplot(data = pwmMat_df_for_ggheatmap, mapping = aes(
    x = positions,
    # Here, 'positions' is the column_name, see previous statement.
    # Do not change it to position_labels
    y = Nucleotides,
    fill = value
  )) +
    ggplot2::geom_tile() +
    ggplot2::theme_bw() +
    ggplot2::xlab(label = "Positions") +
    ggplot2::scale_fill_gradient2(
      name = "", # "Loading",
      low = "white", # "#FFFFFF",
      mid = "white", # FFFFFF",
      high = "#012345"
    ) +
    ggplot2::coord_fixed(ratio = 2.0, clip = "on") +
    ggplot2::theme(
      legend.position = "top",
      legend.justification = "center",
      axis.text.x = element_text(size=rel(0.5), angle = 90, hjust = 1)
    )
  # theme_update(legend.position = "top",
  #              legend.justification = "center"
  #              )

  if (!is.null(savePDFfilename)) {
    if (file.exists(savePDFfilename)) {
      warning("File exists, will overwrite")
    }
    ggplot2::ggsave(p1, device = "pdf", width = 20, height = 2.5)
  }
  return(p1)
}
# test
# position_labels <- seq(5)
# pwmMat <- matrix(rnorm(20), nrow=4)
# rownames(pwmMat) <- c('a', 'c', 'g', 't')
# pwmMat_df <- as.data.frame(pwmMat)
# #
# pwmMat_df <- mutate(pwmMat_df, Nucleotides = rownames(pwmMat_df))
# colnames(pwmMat_df) <- c(position_labels, "Nucleotides")
# #
# pwmMat_df_for_ggheatmap <- melt(pwmMat_df, id.vars=c("Nucleotides"),
# variable.name = "positions")
#
# p1 <- ggplot(data = pwmMat_df_for_ggheatmap, mapping = aes(x = positions,
#                                        y = Nucleotides,
#                                        fill = value)) +
#                       geom_tile() +
#                       theme_bw() +
#                       xlab(label = "Positions")
# print(p1)
#
