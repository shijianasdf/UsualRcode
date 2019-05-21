###################################
#'boxplot单组或者多组数据的分布图
###################################

box.jitter.ggpubr.plot <- function(dat,filepath){
  # @param dat:可视化数据，Var2和value两个字段
  suppressMessages(library(ggpubr))
  suppressMessages(library(ggsignif))
  pdf(filepath)
  p <- ggboxplot(dat, x="Var2", y="value", color="Var2", 
            add = "jitter") +
            rotate_x_text(angle = 45)+
            geom_hline(yintercept = mean(dat$value), linetype = 2)+ # Add horizontal line at base mean
            #stat_compare_means(label = "p.format", method = "wilcox.test")+
            geom_signif(comparisons = list(unique(as.character(dat$Var2))), test = wilcox.test, step_increase = 0.2)
  print(p)
  dev.off()
  return(p)
}

box.jitter.plot <- function(dat,filepath){
  # @param dat:可视化数据，Var2和value两个字段
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggsignif))
  pdf(filepath)
  p <- ggplot(data=dat,mapping = aes(x=Var2,y=value,colour=Var2))+
          geom_boxplot()+
          geom_jitter(shape=16, size=1,position=position_jitter(0.2))+
          labs(x="clonality",y="score value",colour="clonality")+
          geom_signif(comparisons = list(unique(as.character(dat$Var2))), test = wilcox.test, step_increase = 0.2)+
          theme_bw()
  print(p)
  dev.off()
  return(p)
}
