library(ggplot2)
library(scales)#引入percent函数
library(RColorBrewer)
col <- brewer.pal(9,"Blues")[c(2,4,6,8)] #表示使用Blues的2,4,6,8颜色,这样可以区分开
pie.df <- data.frame(group = c("unmutated", "clonal-subclonal", "clonal","subclonal"),
                     value = c(536-513, 283, 224,6))
pie<- ggplot(pie.df, aes(x="", y=value, fill=factor(group,levels = c("unmutated", "subclonal","clonal","clonal-subclonal"))))+ #fill设置因子调颜色顺序
  geom_bar(width = 1, stat = "identity")+
  coord_polar(theta ="y", start=0)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )      
pie.out <- pie + scale_fill_brewer("Blues") + blank_theme +
  theme(axis.text.x=element_blank())+
  labs(title = "% cases based on mutations:",x = NULL, y = NULL, fill = NULL)+
  geom_text(aes(y = value/4 + c(0, cumsum(value)[-length(value)]), 
                label = percent(value/sum(value))),size=5) #调文本位置,暂时没办法，用AI调吧
ggsave(pie.out,filename = "D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.DriverGeneClone/pie.pdf")
