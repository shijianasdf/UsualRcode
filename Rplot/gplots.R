library(gplots)
heatmap.2 (x,

           # dendrogram control
           Rowv = TRUE, 
           Colv=if(symm)"Rowv" else TRUE,
           distfun = dist,
           hclustfun = hclust,
           dendrogram = c("both","row","column","none"),
           reorderfun = function(d, w) reorder(d, w),
           symm = FALSE,

           # data scaling
           scale = c("none","row", "column"),
           na.rm=TRUE,

           # image plot
           revC = identical(Colv, "Rowv"),
           add.expr,

           # mapping data to colors
           breaks, #数值向量，映射数据到颜色,长度比颜色多1
           symbreaks=any(x < 0, na.rm=TRUE) || scale!="none",

           # colors
           col="heat.colors", #颜色向量，标记矩阵颜色，比breaks长度小1。

           # block sepration
           colsep, #数值向量，表记按哪些行分开
           rowsep,
           sepcolor="white",
           sepwidth=c(0.05,0.05),

           # cell labeling
           cellnote,  #数值矩阵，展示矩阵数值
           notecex=1.0,
           notecol="cyan",
           na.color=par("bg"),

           # level trace
           trace=c("column","row","both","none"), #一般设置为none，把那条线隐藏掉
           tracecol="cyan",
           hline=median(breaks),
           vline=median(breaks),
           linecol=tracecol,

           # Row/Column Labeling
           margins = c(5, 5),
           ColSideColors, #和矩阵列数相同长度的向量，储存颜色值，用于标记label
           RowSideColors, #和矩阵行数相同长度的向量，储存颜色值，用于标记label
           cexRow = 0.2 + 1/log10(nr),
           cexCol = 0.2 + 1/log10(nc),
           labRow = NULL, #向量，行名
           labCol = NULL, #向量，列名
           srtRow = NULL,
           srtCol = NULL,
           adjRow = c(0,NA),
           adjCol = c(NA,0),
           offsetRow = 0.5,
           offsetCol = 0.5,
           colRow = NULL,
           colCol = NULL,

           # color key + density info
           key = TRUE, #是否展示legend
           keysize = 1.5, #legend大小
           density.info=c("histogram","density","none"),
           denscol=tracecol,
           symkey = any(x < 0, na.rm=TRUE) || symbreaks,
           densadj = 0.25,
           key.title = NULL,
           key.xlab = NULL,
           key.ylab = NULL,
           key.xtickfun = NULL,
           key.ytickfun = NULL,
           key.par=list(),

           # plot labels
           main = NULL,
           xlab = NULL,
           ylab = NULL,

           # plot layout
           lmat = NULL,
           lhei = NULL,
           lwid = NULL,

           # extras
           extrafun=NULL,
           ...
           )
		  heatmap.2(x,trace="none",key=F,breaks=seq(min(as.numeric(x)), max(as.numeric(x)), length.out=20),col=colorpanel(n=19, low="green", high="red"),ColSideColors=colors1,RowSideColors=colors2)
		   