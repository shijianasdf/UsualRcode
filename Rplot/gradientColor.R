## 取渐变色
gradientColors <- function(colors, resolution=100) {
	#colors, 输入渐变的颜色范围，长度至少为2的向量
	#resolution：数值越大颜色越精细
	#产生主色渐变颜色函数
	ccc <- colorRampPalette(colors, space='rgb', bias=1) # bias>1 渐变色倾向于后面的颜色
	cols <- ccc(resolution)
	return(cols)
}