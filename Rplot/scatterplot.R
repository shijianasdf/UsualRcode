library(ggplot2)
logFC_cutoff <- with(DEG, mean(abs(logFC ) ) + 2 * sd( abs( logFC ) ) )
logFC_cutoff
logFC_cutoff = 1 ## 文章中设置为1
this_tile <- paste0( 'Cutoff for logFC is ', round( logFC_cutoff, 3 ),                       '
                       The number of up gene is ', nrow(DEG[DEG$lab =='up', ] ),                      '
                       The number of down gene is ', nrow(DEG[DEG$lab =='down', ] ) )
volcano <- ggplot(data = DEG, mapping=aes( x = logFC, y = -log10(P.Value), color = lab)) +
                    geom_point(alpha = 0.4, size = 1.75,mapping = aes(colour = lab)) +
                    theme_set( theme_set( theme_bw( base_size = 15 ) ) ) +
                    xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
                    ggtitle( this_tile ) + theme( plot.title = element_text( size = 15, hjust = 0.5)) +
                    scale_colour_manual(values = c('blue','black','red') )
print(volcano) 


volcano.Deseq2.plot <- function(dat.plot,gene,filepath){
  suppressMessages(library(ggplot2))
  pdf(filepath,width=10) #"F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/4.刻画克隆性相关基因/CASP8.volcano.pdf"
  this_tile <- paste0(gene,                    
                      'The number of up gene is ', nrow(dat.plot[dat.plot$lab =='UP',]),                      
                      'The number of down gene is ', nrow(dat.plot[dat.plot$lab =='DOWN',]))
  # Make a basic ggplot2 object with x-y values
  vol <- ggplot(dat.plot, aes(x = log2FoldChange, y = -log10(padj), color = lab))+  
              ggtitle(label = this_tile, subtitle = "Colored by fold-change direction") +
              geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
              #geom_text()
              scale_color_manual(name = "Directionality",
                                 values = c(UP = "#008B00", DOWN = "#CD4F39", NO = "darkgray")) +
              theme_bw(base_size = 14) + # change overall theme
              theme(legend.position = "right") + # change the legend
              xlab(expression(log[2]("FoldChange"))) + # Change X-Axis label
              ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
              geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
              scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
  print(vol)
  dev.off()
  return(vol)
}
