library(pheatmap)
pheatmap(mat, color = colorRampPalette(rev(brewer.pal(n = 7, name =
  "RdYlBu")))(100), kmeans_k = NA, breaks = NA, border_color = "grey60",
  cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE,
  cluster_cols = TRUE, clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean", clustering_method = "complete",
  clustering_callback = identity2, cutree_rows = NA, cutree_cols = NA,
  treeheight_row = ifelse((class(cluster_rows) == "hclust") || cluster_rows,
  50, 0), treeheight_col = ifelse((class(cluster_cols) == "hclust") ||
  cluster_cols, 50, 0), legend = TRUE, legend_breaks = NA,
  legend_labels = NA, annotation_row = NA, annotation_col = NA,
  annotation = NA, annotation_colors = NA, annotation_legend = TRUE,
  annotation_names_row = TRUE, annotation_names_col = TRUE,
  drop_levels = TRUE, show_rownames = T, show_colnames = T, main = NA,
  fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize,
  display_numbers = F, number_format = "%.2f", number_color = "grey30",
  fontsize_number = 0.8 * fontsize, gaps_row = NULL, gaps_col = NULL,
  labels_row = NULL, labels_col = NULL, filename = NA, width = NA,
  height = NA, silent = FALSE, ...);
 #mat 数值矩阵 
 #breaks c(-2,-1,0,1,2) 长度为5，要比color向量长1，将矩阵中的值和颜色映射上
 #color  colorRampPalette(c("white","blue"))(4) 长度为4，设置颜色向量
 #kmeans_k
 #border_color 设置边框的颜色 NA
 #cellwidth 设置每个小格的宽度
 #cluster_rows 是否对行进行聚类
 #clustering_distance_rows 对行聚类的测度
 #cutree_rows 按行分成几个部分
 #treeheight_row 行聚类的树图高度，如果为0，则聚类不显示系统树
 #legend  是否设置legend
 #legend_breaks 向量，对legend设置数值
 #legend_labels 向量，修改legend标签
 #annotation_row  数据框，对行加入其它注释信息
 #display_numbers 是否展示矩阵中的值,也可以自己设置矩阵
 #number_format 矩阵值展示小数点几位
 #number_color 矩阵值展示颜色
 #fontsize 字体大小
library(pheatmap);
tt <- t(cbind(c(1,2,0,0,1,0,0,1,0,1,2,1,0,0,2,1,0),c(2,1,0,0,1,0,0,1,0,1,2,1,0,0,1,1,0)));  
pheatmap(tt,color=c("gray","red","green"),breaks=c(-1,0,1,2),cluster_rows = F,cluster_cols = F,cellwidth = 5, cellheight = 15,border_color = "white");
#(-1,0] (0,1] (1,2] 
  
  
  
  