

#' Visualize a given (PWM) matrix as a sequence logo
#'
#' @param pwmMat Matrix (usually a PWM, but can be any non-normalized matrix) to
#'  be represented as a sequence logo.
#' @param plotMethod For \code{ggseqlogo}; either of "custom", "bits", or "probability".
#' Default is "custom".
#' @param position_labels Labels of the positions in the sequences.
#' @param savePDFfilename Name of the file which will be saved as PDF.
#'
#' @return A ggplot object so you can simply call \code{print} or \code{save}
#' on it later. If \code{savePDFfilename} is given, it is also saved and the
#' \code{ggplot} object returned.
#' @export
#'
#' @import ggplot2
#' @import ggseqlogo
plot_ggseqlogo <- function(pwmMat,
                           plotMethod = "custom",
                           position_labels = NULL,
                           savePDFfilename = NULL) {
  # require(ggplot2)
  # require(ggseqlogo)
  #
  if (!is.matrix(pwmMat)) {
    stop("Expecting a matrix with 4 rows")
  }
  if (sum(dim(pwmMat)) == 2 && is.na(pwmMat)) {
    stop("Empty matrix")
  }
  if (!(nrow(pwmMat) == 4)) {
    stop("Expecting a matrix with 4 rows corresponding to DNA alphabet")
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
  p1 <- ggplot() + geom_logo(pwmMat, method = plotMethod, seq_type = "dna") +
    theme_logo() +
    # theme_update(plot.margin=grid::unit(c(0,0,0,0), "mm"),
    #       axis.ticks.length = unit(0, "lines"),
    #       axis.line = element_blank(),
    #       axis.title.y = element_blank(),
    #       axis.text.y = element_blank()
    #       ) +
    # theme_bw() +
    scale_x_continuous(
      breaks = 1:ncol(pwmMat), labels = position_labels,
      expand = expand_scale(mult = c(0, 0))
    ) +
    theme(axis.text.x = element_text(size=rel(0.5), angle = 90, hjust = 1)) 
    # +
    # ylim(0.0, 2.0)
  #
  if (!is.null(savePDFfilename)) {
    if (file.exists(savePDFfilename)) {
      warning("File exists, will overwrite")
    }
    ggsave(p1, device = "pdf", width = 20, height = 2.5)
  }
  return(p1)
}
