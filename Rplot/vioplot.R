library("vioplot")
vioplot(as.vector(t(na.omit(conservationScore[[4]][6]))),as.vector(t(na.omit(conservationScore[[5]][6]))),as.vector(t(na.omit(conservationScore[[6]][6]))),na.omit(randomregion_phastCons[[2]]),
         names=c("common","non","specific","random"),col=c("red","yellow","blue"));
		 
		 
library(tidyverse)
#devtools::install_github(repo = "IndrajeetPatil/ggstatsplot")
library(ggstatsplot)
library(ggplot2)
library(dplyr)

#函数”%||%”如果a为null返回b，否则返回a
"%||%" <- function(a, b) {
  if (!is.null(a))
    a
  else
    b
}

#利用ggplot2中的ggproto函数自定义几何图形函数GeomFlatViolin，并且让GeomFlatViolin继承自Geom
GeomFlatViolin <-
  ggproto(
    "GeomFlatViolin",  #为函数定义类名，也可以不设置
     Geom,   # GeomFlatViolin继承自Geom
     setup_data = function(data, params) {
      data$width <- data$width %||%
        params$width %||% (resolution(data$x, FALSE) * 0.9)
      
# ymin, ymax, xmin, and xmax define the bounding rectangle for each group
      data %>%
        dplyr::group_by(.data = ., group) %>%
        dplyr::mutate(
          .data = .,
          ymin = min(y),
          ymax = max(y),
          xmin = x,
          xmax = x + width / 2
        )
    },
    
    draw_group = function(data, panel_scales, coord)
    {
      # Find the points for the line to go all the way around
      data <- base::transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
      
      # Make sure it's sorted properly to draw the outline
      newdata <-
        base::rbind(
          dplyr::arrange(.data = base::transform(data, x = xminv), y),
          dplyr::arrange(.data = base::transform(data, x = xmaxv), -y)
        )
      
      # Close the polygon: set first and last point the same
      # Needed for coord_polar and such
      newdata <- rbind(newdata, newdata[1,])
      
      ggplot2:::ggname("geom_flat_violin",
                       GeomPolygon$draw_panel(newdata, panel_scales, coord))
    },
    
    draw_key = draw_key_polygon,
    
    default_aes = ggplot2::aes(
      weight = 1,
      colour = "grey20",
      fill = "white",
      size = 0.5,
      alpha = NA,
      linetype = "solid"
    ),
    
    required_aes = c("x", "y")
  )

#最终定义的函数geom_flat_violin，以后就是调用该函数来完成绘图
geom_flat_violin <-
  function(mapping = NULL, #“映射”则确定如何使用这些数据
           data = NULL,  #可视化的数据，一般为数据框
           stat = "ydensity", #对数据使用的统计学转换
           position = "dodge", #对位置的调整
           trim = TRUE,  #逻辑值 小提琴图是否削掉尾巴
           scale = "area", #让所有小提琴有相同的面积           show.legend = NA,   #是否展示图例
           inherit.aes = TRUE,  #默认T，继承映射（mapping）
           ...) {
    ggplot2::layer( #ggplot2 layer函数定义一个图层(一般用geom_*定义)
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomFlatViolin, #上面定义的函数，用自定义函数来可视化数据
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(trim = trim,
                    scale = scale,
                    ...)
    )
  }


#调用上面定义函数geom_flat_violin，完成小提琴图
ggplot2::ggplot(data = iris, mapping = ggplot2::aes(x = Species, y = Sepal.Length)) + geom_flat_violin() 
#@param  data: 一个数据框,要可视化的数据。
#@param  mapping: ”映射”,确定如何使用这个数据框,比如这里将iris中的Species字段映射到X轴(有3种值setosa,versicolor,virginica)，Sepal.Length字段映射到Y轴。
#@param  geom_flat_violin: 利用上面自定义小提琴图型函数绘制新的图层，告诉R我们要画小提琴图，具体参数见geom_flat_violin。


ggplot2::ggplot(data = iris, mapping = ggplot2::aes(x = Species, y = Sepal.Length, fill = Species)) + geom_flat_violin(trim = FALSE) 
#@param  data: 一个数据框,要可视化的数据。
#@param  mapping: ”映射”,确定如何使用这个数据框,比如这里将iris中的Species字段映射到X轴(有3种值setosa,versicolor,virginica)，Sepal.Length字段映射到Y轴，Species映射到fill填充色上。
#@param  geom_flat_violin:利用上面自定义小提琴图型函数绘制新的图层，告诉R我们要画小提琴图,作图数据没有设置，默认来自于iris，映射mapping也没有设置，默认同ggplot mapping参数相同，在该函数中trim=FALSE表示留下小提琴的尾巴，具体参数见geom_flat_violin。


# 对生成的小提琴图继续加入新的图层，所有新的图层如果没有指定data和mapping参数，那么默认data和mapping参数继承于ggplot函数的data和mapping
ggplot2::ggplot(data = iris,
                mapping = ggplot2::aes(x = Species, y = Sepal.Length, fill = Species)) +
#自定义小提琴图层 
  geom_flat_violin(scale = "count", trim = FALSE，show.legend=T) +
# stat_summary加入统计信息图层，比如中位数等，具体参数?stat_summary
  ggplot2::stat_summary(    
fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "pointrange",
    position = ggplot2::position_nudge(x = 0.05, y = 0)
  ) +
#加入点图图层  具体参数?geom_dotplot
  ggplot2::geom_dotplot(
    binaxis = "y",
    dotsize = 0.5,
    stackdir = "down",
    binwidth = 0.1,
    position = ggplot2::position_nudge(x = -0.025, y = 0)
  ) +
#加入坐标轴图层 具体参数?labs
  ggplot2::labs(x = "Species", y = "Sepal length (cm)") + 
#定义主题(我理解为样式，整体风格，背景色等)来自ggstatsplot
  ggstatsplot::theme_mprl() + 
#定义主题(我理解为样式，整体风格，背景色等)来自ggplot2
  ggplot2::theme(legend.position = "none")
