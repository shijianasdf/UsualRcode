##stack barplot
head(df,9)
# gene            group     cases
# ACVR1B.1 ACVR1B           Clonal 77.777778
# ACVR1B.2 ACVR1B        Subclonal 18.518519
# ACVR1B.3 ACVR1B Clonal-Subclonal  3.703704
# AKT1.1     AKT1           Clonal 91.666667
# AKT1.2     AKT1        Subclonal  8.333333
# AKT1.3     AKT1 Clonal-Subclonal  0.000000
# ANK1.1     ANK1           Clonal 78.048780
# ANK1.2     ANK1        Subclonal 19.512195
# ANK1.3     ANK1 Clonal-Subclonal  2.439024
library(ggplot2)
library(RColorBrewer)
col <- brewer.pal(9,"Blues")[2,4,6,8] #表示使用Blues的2,4,6,8颜色
stack.barplot <- ggplot(data=df, mapping=aes(x=gene, y=cases, fill=factor(group,levels=c("Subclonal","Clonal","Clonal-Subclonal"))))+ #levels=c("Subclonal","Clonal","Clonal-Subclonal" 调整颜色顺序
                        geom_bar(stat="identity",width=0.7,position="stack")+ #default width 0.9 position="stack" position=position_dodge(0.8)
                        scale_fill_manual(values = col[2:4])+
                        theme_minimal()+ labs(fill="clone status",y="% cases")+
                        #x坐标文本字体角度angle，vjust垂直移动多少，hjust水平移动多少，size设置字体大小正常的0.6倍，调整图中panel的网格为空
                        theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.5,size = rel(0.4)),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
                        theme(legend.position = "top",legend.key.size = unit(0.5, "cm"))+  #调整图例位置，大小
                        theme(legend.title=element_text(size=9),legend.text=element_text(size=9))+  #调整图例标题大小 调整图例内部字体大小
ggsave(stack.barplot,filename = "D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/1.DriverGeneClone/stack.barplot.pdf")

##group barplot
library(ggplot2)
picture <- ggplot(data=picture.three.data,mapping=aes(genes,number,fill=status))+ #数据框 
                xlim(0, 80)+
                geom_bar(stat="identity", width = 0.5, position = position_dodge(0.9)) + #柱的宽度
                scale_fill_manual(values = c("RoyalBlue", "red"), labels=c("Clonal", "Subclonal"))+
                geom_text(mapping = aes(label= number))+ #加文本
                theme_set(theme_bw())+ #背景白色
                theme(panel.grid.major=element_line(colour=NA))+ #画板panel的格线去掉
                theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 90))+
                coord_flip() #坐标翻转
