#' Plot a heatmap using ggplot2.
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
#' @examples
#'
#'
plot_ggheatmap <- function(pwmMat, position_labels=NULL, savePDFfilename=NULL){
        # require(ggplot2)
        # require(reshape2)
        pwmMat_df <- as.data.frame(pwmMat)
        #
        pwmMat_df <- mutate(pwmMat_df, Nucleotides = rownames(pwmMat_df))
        colnames(pwmMat_df) <- c(position_labels, "Nucleotides")
        #
        pwmMat_df_for_ggheatmap <- melt(pwmMat_df, id.vars=c("Nucleotides"), variable.name = "positions")
        p1 <- ggplot(data = pwmMat_df_for_ggheatmap, mapping = aes(x = positions,
                                                                   y = Nucleotides,
                                                                   fill = value)) +
                    geom_tile() +
                    theme_bw() +
                    xlab(label = "Positions") +
                    scale_fill_gradient2(name = "", #"Loading",
                                         low = "white", #"#FFFFFF",
                                         mid = "white", #FFFFFF",
                                         high = "#012345"
                                        ) +
                    theme(legend.position = "top",
                          legend.justification = "center"
                    )
                    # theme_update(legend.position = "top",
                    #              legend.justification = "center"
                    #              )

        if(!is.null(savePDFfilename)){
                ggsave(p1, device="pdf", width=20, height=2.5)
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
  # pwmMat_df_for_ggheatmap <- melt(pwmMat_df, id.vars=c("Nucleotides"), variable.name = "positions")
  #
  # p1 <- ggplot(data = pwmMat_df_for_ggheatmap, mapping = aes(x = positions,
  #                                        y = Nucleotides,
  #                                        fill = value)) +
  #                       geom_tile() +
  #                       theme_bw() +
  #                       xlab(label = "Positions")
  # print(p1)
  #


