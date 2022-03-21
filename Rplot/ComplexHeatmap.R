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

#重新绘图，加入样本信息
    infercnv.reference.matrix <- fread("/Users/biofly/project/shijian/Results/inferCNV/6_sample_inferCNV/infercnv.references.txt")
    infercnv.observations.matrix <- fread("/Users/biofly/project/shijian/Results/inferCNV/6_sample_inferCNV/infercnv.observations.txt")
    infercnv.reference.matrix <- infercnv.reference.matrix %>% as.data.frame() %>% tibble::column_to_rownames("V1") 
    infercnv.observations.matrix <- infercnv.observations.matrix %>% as.data.frame() %>% tibble::column_to_rownames("V1") 
    match(rownames(sce_merge_res$sce@meta.data),colnames(infercnv.reference.matrix))
    #正常细胞heatmap
    reference.metadata <- sce_merge_res$sce@meta.data[match(colnames(infercnv.reference.matrix),rownames(sce_merge_res$sce@meta.data)),]
    col_celltype <- palettes(category = "box", #创建颜色画板
             palette = "paired1",
             alpha = 1,
             counts = 50, 
             show_col = TRUE, 
             show_message = FALSE)
    col_celltype1 <- palettes(category = "box", 
                             palette = "nrc",
                             alpha = 1,
                             counts = 50, 
                             show_col = TRUE, 
                             show_message = FALSE)
    col_orig.ident <- palettes(category = "box", 
                              palette = "jama",
                              alpha = 1,
                              counts = 50, 
                              show_col = TRUE, 
                              show_message = FALSE)
    normal_cell_anno <- rowAnnotation(orig.ident = reference.metadata$orig.ident,
                                      celltype=reference.metadata$celltypes,
                                      col=list(orig.ident = c("P1_653303"=col_orig.ident[1],"P2_656098"=col_orig.ident[2],"P3_651038"=col_orig.ident[3],"P4_662693"=col_orig.ident[4],"P5_665746"=col_orig.ident[5],"P6_656098"=col_orig.ident[6]),
                                               celltype = c("immue cell_25"=col_celltype[1],"immue cell_23"=col_celltype[2],  "endothelial_17"=col_celltype[3],
                                                            "immue cell_15"=col_celltype[4],  "fibroblasts_13"=col_celltype[5], "immue cell_10"=col_celltype[6],  "immue cell_8"=col_celltype[7],
                                                            "immue cell_7"=col_celltype[8],   "immue cell_6"=col_celltype1[1],   "immue cell_5"=col_celltype1[2],   "immue cell_1"=col_celltype1[3],
                                                            "immue cell_0"=col_celltype1[4])))
    cols <- colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(c(-5,0,5)))
    col <- circlize::colorRamp2(c(-5,0,5), cols)
    p <- Heatmap(t(t(scale(t(infercnv.reference.matrix)))), #可视化矩阵（行名以及列名）
                 show_column_names = F, 
                 show_row_names = F,
                 border = T, #
                 row_split = reference.metadata$celltypes, #一个向量，基于哪个指标对矩阵进行拆分，从而基于拆分的子矩阵进行聚类
                 row_dend_width = unit(20, "mm"),
                 row_title = "cells",
                 col=col, #设置颜色 由颜色映射函数circlize::colorRamp2(c(-5,0,5), cols)生成
                 cluster_columns = F,
                 left_annotation=normal_cell_anno)  #左侧注释bar 由rowAnnotation(orig.ident = reference.metadata$orig.ident,celltype=reference.metadata$celltypes)生成
    draw(p)
    #肿瘤细胞heatmap
    observations.metadata <- sce_merge_res$sce@meta.data[match(colnames(infercnv.observations.matrix),rownames(sce_merge_res$sce@meta.data)),]
    tumour_cell_anno <- rowAnnotation(orig.ident = observations.metadata$orig.ident,
                                      celltype=observations.metadata$celltypes,
                                      col=list(orig.ident = c("P1_653303"=col_orig.ident[1],"P2_656098"=col_orig.ident[2],"P3_651038"=col_orig.ident[3],"P4_662693"=col_orig.ident[4],"P5_665746"=col_orig.ident[5]),
                                               celltype = c("EPCAM cell_2"=col_celltype[1],"EPCAM cell_12"=col_celltype[2],  "EPCAM cell_4"=col_celltype[3],
                                                            "hepatocytes_24"=col_celltype[4],  "EPCAM cell_3"=col_celltype[5], "EPCAM cell_18"=col_celltype[6],  "EPCAM cell_20"=col_celltype[7],
                                                            "EPCAM cell_19"=col_celltype[8],   "EPCAM cell_9"=col_celltype1[1],   "EPCAM cell_11"=col_celltype1[2],   "EPCAM cell_14"=col_celltype1[3],
                                                            "EPCAM cell_22"=col_celltype1[4], "EPCAM cell_16"=col_celltype1[5],"EPCAM cell_21"=col_celltype1[6])))
    # col_scale <- rev(RColorBrewer::brewer.pal(9, "RdBu"))
    # col_palette <- grDevices::colorRampPalette(col_scale)(n = 50)
    cols <- colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(c(-3,0,3)))
    col <- circlize::colorRamp2(c(-3,0,3), cols)
    p <- Heatmap(t(t(scale(t(infercnv.observations.matrix)))),
                 show_column_names = F, 
                 show_row_names = F,
                 border = T,
                 row_split = observations.metadata$celltypes,
                 row_dend_width = unit(20, "mm"),
                 row_title = "cells",
                 col=col,
                 cluster_columns = F,
                 left_annotation=tumour_cell_anno) 
    draw(p)
