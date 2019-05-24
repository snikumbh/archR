viz_all_factors_in_combined_heatmaps_seqlogos <- function(featuresMatrix, position_labels=NA, add_pseudo_counts = F, savePDFfilename=NULL){
  suppressMessages( require(cowplot) )
  suppressMessages( require(gridExtra) )
  #
  invisible(apply(featuresMatrix, MARGIN = 2, function(x){
                      pwm <- make_sinuc_PWMs(x, add_pseudo_counts = add_pseudo_counts, scale=F)
                      #
                      # Heatmap on top
                      p1 <- plot_ggheatmap(pwmMat=pwm, 
                                           position_labels=positions, 
                                           savePDFfilename=savePDFfilename
                      )
                      # Make adjustments for alignment
                      # p1 + theme(legend.position = "top", 
                      #                   legend.justification = "center"
                      # )
                      # Seqlogo below
                      p2 <- plot_ggseqlogo(pwmMat=pwm, 
                                           position_labels=positions, 
                                           savePDFfilename=savePDFfilename
                      )
                      # Make adjustments for alignment
                      p2 + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")#,
                                        #axis.ticks.length = unit(0.1, "lines"),
                                        #axis.line.y = element_blank(),
                                        #axis.title.y = element_blank(),
                                        #axis.text.y = element_blank()
                      )
                      
                      #do.call(gridExtra::grid.arrange, c(p_list, ncol=2))
                      gridExtra::grid.arrange(p1, p2)
                      #print(plot_grid(p1, p2,  ncol = 1, align = 'h'))
                    }
          ))
}