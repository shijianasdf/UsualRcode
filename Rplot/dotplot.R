vignette("ggplot2-specs")
head(colorectal.select.table)
# patient                   mutation_id Hugo_Symbol
# 25  TCGA-3L-AA1B TCGA-3L-AA1B-01:14:65544658:G         MAX
# 62  TCGA-3L-AA1B TCGA-3L-AA1B-01:1:115258747:C        NRAS
# 111 TCGA-3L-AA1B TCGA-3L-AA1B-01:3:178936091:G      PIK3CA
# 126 TCGA-3L-AA1B TCGA-3L-AA1B-01:5:112173917:C         APC
# 127 TCGA-3L-AA1B TCGA-3L-AA1B-01:5:112175639:C         APC
# 147 TCGA-3L-AA1B TCGA-3L-AA1B-01:7:140453146:G        BRAF
# Variant_Classification absolute.ccf CI95.timing
# 25       Missense_Mutation         1.00      Clonal
# 62       Missense_Mutation         1.00      Clonal
# 111      Missense_Mutation         1.00      Clonal
# 126      Nonsense_Mutation         1.00      Clonal
# 127      Nonsense_Mutation         0.56   Subclonal
# 147      Missense_Mutation         1.00      Clonal
# prob.clonal.timing
# 25              Clonal
# 62              Clonal
# 111             Clonal
# 126             Clonal
# 127          Subclonal
# 147             Clonal
library(ggplot2)
library(RColorBrewer)
col <- brewer.pal(9,"Blues")[2,4,6,8] #表示使用Blues的2,4,6,8颜色
dot.plot<-ggplot(data=colorectal.select.table, mapping=aes(x=Hugo_Symbol, y=absolute.ccf,fill=factor(CI95.timing,levels = c("Subclonal","Clonal")))) + 
                  geom_jitter(shape= 21,alpha = 1,size=1,width=0.2,color="black",stroke=0.4)+ #You'll have to use shapes from 21 to 25. These are the ones that have colour and fill properties
                  labs(fill = "Legend Title")+ #设置legend题目为"Legend Title"
                  theme_bw()+scale_fill_manual(values=col[2:3],labels = c("Subclonal","Clonal")) +  #要先设置默认的主题风格
                  theme(legend.position = "top",axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.5,size = rel(0.5)),panel.grid.major = element_line(size = 0.1),panel.grid.minor = element_line(size = 0.1))
dot.plot
ggsave(dot.plot,filename = "D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.DriverGeneClone/dotplot.pdf")

dot.plot <- ggplot(data=colorectal.select.table, mapping=aes(x=Hugo_Symbol, y=as.numeric(absolute.ccf),fill=factor(CI95.timing,levels = c("Subclonal","Clonal")))) + 
                  geom_dotplot(shape=21,binaxis = "y",dotsize = 0.4,width=0.05,position = "jitter",color="black")+theme_bw()+scale_fill_manual(values = col[2:3])+
                  theme(legend.position = "top",axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.5,size = rel(0.5)))
dot.plot


