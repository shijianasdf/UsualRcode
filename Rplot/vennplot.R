library(VennDiagram)
venn.plot <- venn.diagram(
        x = list(
                diff.genes = c(sig.down.genes.FDR.0.05.symbols,sig.up.genes.FDR.0.05.symbols),
                enhancer.target.genes =breast_specific_enhancers_target_genes
                ),
        filename = NULL,
		fill = c("dodgerblue", "goldenrod1"),
		cat.col = c("dodgerblue", "goldenrod1")
        );
grid.draw(venn.plot);








venn.diagram(x, filename, height = 3000, width = 3000, resolution = 500,units = "px", 
				compression = "lzw", na = "stop", main = "", sub = "",main.pos = c(0.5, 1.05), 
				main.fontface = "plain", main.fontfamily = "serif",main.col = "black", main.cex = 1, 
				main.just = c(0.5, 1),sub.pos = c(0.5, 1.05), sub.fontface = "plain", sub.fontfamily = "serif",
				sub.col = "black", sub.cex = 1, sub.just = c(0.5, 1),category.names = names(x), 
				force.unique = TRUE, ...)
#x:一个向量列表
#fill:韦恩图中每个圆填充的颜色
#col:韦恩图中每个圆边缘颜色的设置
#label.col:韦恩图中数字的颜色
#lwd:整数型参数，韦恩图中每个圆边线的粗细
#lty:韦恩图中每个圆边缘线的类型，例如lty="dotted",边缘线为虚线，默认情况下是实线
#cex:整数型参数，韦恩图中数字的大小
#main.cex与sub.cex:韦恩图中主副标题字的大小
#main:字符型参数，韦恩图的主标题
#sub:字符型参数，韦恩图的副标题
#ratation.degree:整数型参数，韦恩图旋转的角度
#cat.cex:标签名字的大小
#cat.col:标签名字的颜色
#cat.pos:标签名字相对于圆的位置
#cat.dist:标签名字相对于圆的角度