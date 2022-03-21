#基因组突变热图
library(ComplexHeatmap)
mat = read.table(paste0(system.file("extdata", package = "ComplexHeatmap"), 
                        "/tcga_lung_adenocarcinoma_provisional_ras_raf_mek_jnk_signalling.txt"), 
                 header = TRUE,stringsAsFactors=FALSE, sep = "\t")
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat=  mat[, -ncol(mat)]
mat = t(as.matrix(mat))
mat[1:3, 1:3]
# TCGA-05-4384-01 TCGA-05-4390-01 TCGA-05-4425-01
# KRAS "  "            "MUT;"          "  "           
# HRAS "  "            "  "            "  "           
# BRAF "  "            "  "            "MUT;AMP;" 
col = c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
  }
)
#In this case, we need to define a function to extract different alteration types and pass the function to get_type argument. 
#The function should return a vector of alteration types.
oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "OncoPrint for TCGA Lung Adenocarcinoma, genes in Ras Raf MEK JNK signalling",
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "HOMDEL", "MUT"), 
                                      labels = c("Amplification", "Deep deletion", "Mutation")))
dim(mat)

library(ComplexHeatmap)
cols <- colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(c(-5,0,5)))
col <- circlize::colorRamp2(c(-5,0,5), cols)
sam_anno <- HeatmapAnnotation(orig.ident = annotation_col$orig.ident)
p <- Heatmap(t(scale(t(pathway.ssgsea.matrix[sig.pathway[1:183],]))),show_column_names = F, column_split = annotation_col$orig.ident,col=col,top_annotation=sam_anno,row_names_max_width = unit(25, "cm"),row_names_gp = gpar(fontsize = 5,fontface="bold")) 
draw(p,annotation_legend_side = "left")

reference.metadata <- sce_merge_res$sce@meta.data[match(colnames(infercnv.reference.matrix),rownames(sce_merge_res$sce@meta.data)),]
    normal_cell_anno <- rowAnnotation(orig.ident = reference.metadata$orig.ident,celltype=reference.metadata$celltypes) 
    col_scale <- rev(RColorBrewer::brewer.pal(9, "RdBu"))
    col_palette <- grDevices::colorRampPalette(col_scale)(n = 50)
    cols <- colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(c(-5,0,5)))
    col <- circlize::colorRamp2(c(-5,0,5), cols)
    p <- Heatmap(t(t(scale(t(infercnv.reference.matrix)))), #可视化矩阵（行名以及列名）
                 show_column_names = F, 
                 show_row_names = F,
                 border = T, #
                 row_split = reference.metadata$celltypes, #基于哪个指标对矩阵进行拆分，从而基于拆分的子矩阵进行聚类
                 row_dend_width = unit(20, "mm"),
                 row_title = "cells",
                 col=col, #设置颜色 由circlize::colorRamp2(c(-5,0,5), cols)生成
                 cluster_columns = F,
                 left_annotation=normal_cell_anno)  #左侧注释bar 由rowAnnotation(orig.ident = reference.metadata$orig.ident,celltype=reference.metadata$celltypes)生成
    draw(p)

