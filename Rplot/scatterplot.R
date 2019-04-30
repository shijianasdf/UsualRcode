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