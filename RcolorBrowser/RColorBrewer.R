################################################################################
#共有三组颜色可供使用：
#    Sequential，按顺序渐变的(范围：3-9)。 - Light colours for low data, dark for high data Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd
#    Diverging，彼此之间差异变化较大的(范围3-11)。 -  Light colours for mid-range data, low and high contrasting dark colours  BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral 
#    Qualitative，这个用于最大程度地显示不同类之间的差别。 - Colours designed to give maximum visual difference between classes 
################################################################################
library(RColorBrewer)
display.brewer.all()  #显示所有可用颜色
display.brewer.all(type = "div")
display.brewer.all(type = "seq")
display.brewer.all(type = "qual")
col=brewer.pal(8,"Set3") #表示使用Set3的8种颜色
newpalette<-colorRampPalette(brewer.pal(9,"Blues"))(10) #colorRampPalette将Blues系列九个颜色进行了延伸，产生了10种渐变色。
#colorRampPalette生成画板渐变色，输入是一个向量颜色



library(RColorBrewer)
pdf("/IJob/J26/lunwen/5.survival/enhancer_box5.pdf")
col=brewer.pal(3,"Set3")
a <- 5:9
b <- 3:9
d <- 4:6
b1 <- sample(a,159,replace=T)
b2 <- sample(b,126,replace=T) 
b3 <- sample(d,208,replace=T)
boxplot(list(b1,b2,b3),col=col,ylab="log2(RPKM+1)")
dev.off()





library(RColorBrewer)
pdf("/IJob/J26/lunwen/5.survival/enhancer_box5.pdf")
col=brewer.pal(5,"Set3")
a <-c(47,5,52,39,28)
names(a) <- c("BRCA","COAD","GBM","LIHC","LUSC")
boxplot(list(b1,b2,b3),col=col,ylab="log2(RPKM+1)")
dev.off()