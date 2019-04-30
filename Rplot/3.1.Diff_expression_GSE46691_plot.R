#' @author dragon 2017.7.25
#' @description 绘制小提琴图

# ggplot2的绘图过程我将其分为三步，
# 数据处理，也就是所谓的统计方法，
# 数据转换，将数据转换成二维指定的坐标系，将映射至颜色，大小，线型等属性数据转换成属性值。
# 调整坐标，图标，背景等，绘图输出。

####绘制NEAT,HNRNPU-AS1的在转移和非转移中的差异表达图
load("/IJob/J20/Data/metastasis/lzl/prostate_cancer/GSE46691/results/met_exp.RData")
met_exp<-as.data.frame(met_exp)
pdf(file="/IJob/J20/Data/metastasis/lzl/prostate_cancer/GSE46691/plots/met_exp_violin.pdf",pointsize=12)
library(ggplot2)
toplncRNA<-colnames(met_exp)[1:9]
toplncRNA<-toplncRNA[c(2,9)]
lab<-sapply(toplncRNA,function(x){
	a<-paste(x,"_",met_exp$Met_label,sep="")
	return(a)
	})
lab<-as.character(lab)
#数据与图形属性之间的映射关系称为mapping，在ggplot2中用aes()，x一般需求为因子型,即使你不设定因子，函数内部好像也会自动将其转为因子型
#但如果使用因子本身的默认排序，绘图时并不能按照自己想要的顺序排列，即定义了有序因子，可以自定义图形的排列顺序
labb<-ordered(lab, levels = c('ENSG00000245532_0','ENSG00000245532_1','ENSG00000188206_0',"ENSG00000188206_1"))

lnc_exp<-as.matrix(met_exp[,c(2,9)])
lnc_exp<-as.numeric(lnc_exp)
pp<-data.frame(exp=lnc_exp,lab=labb)

p <- ggplot(pp, aes(x=lab,y=exp,fill=lab))
p + geom_violin(trim=F,size=0,show.legend=F)+geom_jitter(shape= 21,alpha = 0.3,size=0.25, show.legend=F)+geom_boxplot(fill="white",width=0.15,size=0.5,outlier.size=0.5)+ggtitle("lncRNA expression")+theme(panel.background = element_rect(fill = "white", colour = "black",size = 1))
dev.off()

