##### Figure plotting #####
## Figure 3A ## 
####Subclusters of Epi-SL
EpiSL<-subset(seurat,idents='Epi-SL')
EpiSL <- FindVariableFeatures(EpiSL, selection.method = "vst", nfeatures = 4000)
top10<-head(VariableFeatures(EpiSL),10)
top10
scale.genes <-  rownames(EpiSL)
EpiSL <- ScaleData(EpiSL, features = scale.genes)
EpiSL<- RunPCA(EpiSL, features = VariableFeatures(object =EpiSL))
plot2 <- ElbowPlot(EpiSL, ndims=20, reduction="pca") 
plot2
pc.num=1:5
EpiSL <- FindNeighbors(EpiSL, dims = pc.num)
EpiSL <- FindClusters(EpiSL, resolution = 0.5)
table(EpiSL@meta.data$seurat_clusters)
metadata <- EpiSL@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)

EpiSL <- RunUMAP(EpiSL, n.neighbors = 50,n.epochs=50,dims = pc.num)
embed_umap <- Embeddings(EpiSL, 'umap')
write.csv(embed_umap,'Epi cluster/EpiSL cluster/embed_umap n.neighbors = 50,n.epochs=50 20211125.csv') 
plot2 = DimPlot(EpiSL,label = T, reduction = "umap") 
plot2
DimPlot(EpiSL,label = T, group.by = "group",reduction = "umap") 
DimPlot(EpiSL,label = T, group.by = "celltype",reduction = "umap") 
ggsave("Epi cluster/EpiSL cluster/UMAP n.neighbors = 50,n.epochs=50 20211125.png", plot = plot2, width = 8, height = 7)
diff.wilcox = FindAllMarkers(EpiSL)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "Epi cell_identify/EpiSL cell_identify/diff_genes_wilcox211125.csv", row.names = F)
write.csv(top10, "Epi cell_identify/EpiSL cell_identify/top10_diff_genes_wilcox211125.csv", row.names = F)

#根据以上注释结果对细胞群添加注释信息
EpiSL@meta.data$celltypes3<-EpiSL@meta.data$seurat_clusters
Idents(EpiSL) <- "celltypes3" 
EpiSL<- RenameIdents(EpiSL, '1'='Epi-1','2'='Epi-1','3'='Epi-1','5'='Epi-1',
                     '7'='Epi-1','0'='Epi-2','4'='Epi-2','6'='Epi-2',
                     '8'='Epi-2','9'='Epi-2')
EpiSL$celltypes3<-Idents(EpiSL)
plot2<-DimPlot(EpiSL, reduction = "umap",cols=c("#33A02C","#E31A1C"),
               group.by ='celltypes3',label = T)
plot2
save(EpiSL,file='seurat_EpiSL.Rdat')
ggsave("Epi cluster/EpiSL cluster/UMAP 20211125.pdf", plot = plot2, width = 6, height = 5)



## Figure 3B ## 
# set colors
blue   <- "#5bc0eb"
yellow <- "#fde74c"
green  <- "#9bc53d"
red    <- "#f25f5c"
purple <- "#531f7a"
grey   <- "#8693ab"
orange <- "#fa7921"
white  <- "#f2d7ee"
darkred   <- "#F2042C"
lightred  <- "#FF7FBF"
lightblue <- "#B2EBFF"
darkblue  <- "#1d00ff"
cherry    <- "#700353"
lightgrey <- "#dcddde"
nake <- "#F8C364"
gold <- "#ECE700"
cyan <- "#00B3D0"
sun  <- "#E53435"
peach  <- "#E43889"
violet <- "#89439B"
soil   <- "#EC7D21"
lightgreen <- "#54B642"
darkblue   <- "#21498D"
darkgreen  <- "#009047"
brown      <- "#874118"
seagreen   <- "#008B8A"
jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")
heatmap.BlWtRd <- c("#6699CC","white","#FF3C38")
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
heatmap.GrWtRd <- c("#2b2d42","#8d99ae","#edf2f4","#ef233c","#d90429")
heatmap.L.BlYlRd <- c("#4281a4","#9cafb7","#ead2ac","#e6b89c","#fe938c")
heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")
heatmap.fancy <- c("#10040A", "#2A0B35", "#4D155B", "#73215B", "#9C3558", "#C34D44", "#E07038", "#F2981C", "#F2CA51", "#FAF6A3")

# customized function
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

gmt2list <- function(annofile){
  if (!file.exists(annofile)) {
    stop("There is no such gmt file.")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are accepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
}


# load R package
library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(tidyr)
library(Seurat)
library(cowplot)
library(ComplexHeatmap)
library(GSVA)
library(MOVICS)
library(gplots)

# load data
load("Epi_metadata.Rdat")
load("Epi_metadata20211125V2.Rdat")
dat <- Epi_metadata[,c("new.patient","group","celltypes4")]
dat$new.patient <- factor(dat$new.patient,levels = rev(c("NC1","NC2","NC3","HP1","HP2","HP3","HP4","SSL1","SSL2","SSL3","SSL4","SSL5","SSLD","TSA1","TSA2","TSA3","TSA4","TSA5")))
dat$celltypes <- factor(dat$celltypes4,levels = rev(c('BEST4+ colonocytes','Mature colonocytes','Goblet cells','Intermediate',
                                                      'Transit amplifying cells','Epi-2','Epi-1','Tuft cells',
                                                      'Enteroendocrine cells')))
# mycol <- rev(viridis(10))
# mycol <- rev(inferno(10))
# mycol <- brewer.pal(10,"Spectral")
mycol <- brewer.pal(10,"Paired")[c(1,3,5,7,9,2,4,6,8,10)]

cell.col <- c("BEST4+ colonocytes" = mycol[1],
              "Mature colonocytes" = mycol[2],
              "Goblet cells" = mycol[3],
              "Intermediate" = mycol[4],
              "Transit amplifying cells" = mycol[5],
              "Epi-2" = "#E31A1C",
              "Epi-1" = "#33A02C",
              #"Epi-3" = mycol[8],
              "Tuft cells" = mycol[9],
              "Enteroendocrine cells" = mycol[10])
save(cell.col,file = "cell.col20211125V2.RData")

dat.bar <- as.data.frame(table(dat$new.patient,dat$celltypes))
dat.bar <- dat.bar%>% group_by(Var1) %>% 
  mutate(Pct = Freq / sum(Freq) * 100) %>% as.data.frame()
colnames(dat.bar) <- c("patient","celltype","count","pct")

# left barplot
p.left <- ggplot(dat.bar, aes(x = patient, y = count,fill = celltype)) +
  scale_fill_manual(values  = cell.col) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(name = "",position = "top") +
  #theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.y = element_line(size = 0.6),
        axis.text.y = element_text(colour = "black",size = 10),
        axis.title.x = element_text(vjust = -0.3,size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.title = element_blank()) +
  coord_flip() +
  scale_y_reverse(expand = c(0.01,0),
                  name = "Epithelial cell number", position = "right")
p.left


p.right <- ggplot(dat.bar, aes(x = patient, y = pct, fill=celltype)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values  = cell.col) +
  scale_x_discrete(name = "") +
  #theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.y = element_line(size = 0.6),
        axis.text.y = element_blank(),
        axis.title.x = element_text(vjust = -0.3,size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.title = element_blank()) +
  coord_flip() +
  scale_y_continuous(expand = c(0.01,0),
                     name = "Epithelial cell proportion (%)", position = "right")

p.right

# middle label
# pp <- ggplot() +
#   geom_text(data = dat.bar,
#              aes(label = patient, x = patient),
#              y = 0.5,
#              color = "black",
#              size = 0.8*11/.pt, # match font size to theme
#              hjust = 0.5, vjust = 0.5) +
#   theme_minimal()+
#   theme(axis.line.y =element_blank(),
#         axis.ticks.y =element_blank(),
#         axis.text.y =element_blank(),
#         axis.title.y =element_blank(),
#         axis.title.x =element_blank(),
#         plot.margin = unit(c(0.3, 0, 0.3, 0), "lines")
#   ) +
#   #guides(fill = FALSE) +
#   coord_flip() +
#   scale_y_reverse()
# pp

# combine figures
# pal <- p.left + pp + p.right +
#   plot_layout(widths = c(7,1,7), guides = 'collect') & theme(legend.position = 'right',legend.key.size = unit(0.5, 'cm'))
# pal

pal <- p.left + p.right +
  plot_layout(widths = c(7,7), guides = 'collect') & theme(legend.position = 'right',legend.key.size = unit(0.5, 'cm'))
pal
ggsave(filename = "distribution of epi cells.pdf", width = 6,height = 3)

dat.bar <- as.data.frame(table(dat$group,dat$celltypes))
dat.bar <- dat.bar%>% group_by(Var1) %>% 
  mutate(Pct = Freq / sum(Freq) * 100) %>% as.data.frame()
colnames(dat.bar) <- c("group","celltype","count","pct")
dat.bar$group <- factor(dat.bar$group,levels = rev(c("TSA","SSLD","SSL","HP","NC")))

p.left <- ggplot(dat.bar, aes(x = group, y = count,fill = celltype)) +
  scale_fill_manual(values  = cell.col) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(name = "",position = "bottom") +
  ylab("Epithelial cell  number") + 
  theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.y = element_line(size = 0.6),
        axis.title = element_text(size = 12,colour = "black"),
        axis.text.y = element_text(colour = "black",size = 10),
        axis.title.x = element_text(vjust = -0.3,size = 12),
        axis.text.x = element_text(size = 10, color = "black",angle = 45,hjust = 1),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.title = element_blank())

p.left

p.right <- ggplot(dat.bar, aes(x = group, y = pct, fill=celltype)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values  = cell.col) +
  scale_x_discrete(name = "") +
  theme_bw() +
  ylab("Epithelial cell proportion (%)") + 
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.y = element_line(size = 0.6),
        axis.title = element_text(size = 12,colour = "black"),
        axis.title.x = element_text(vjust = -0.3,size = 12),
        axis.text.x = element_text(size = 10, color = "black",angle = 45,hjust = 1),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.title = element_blank())
p.right

pal <- p.left + p.right +
  plot_layout(widths = c(6,6), guides = 'collect') & theme(legend.position = 'right',legend.key.size = unit(0.4, 'cm'))
pal
ggsave(filename = "distribution of epi cells2.20211125V2.pdf", width = 5,height = 3)



## Figure 3C ## 
library(monocle)
Idents(seurat) <- "celltypes4" 
serrated_epi<-subset(seurat,idents=c('Transit amplifying cells','Epi-1','Epi-2'))
dir.create("pseudotime")
scRNAsub <- serrated_epi  #scRNAsub是 子集seurat对象
data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)

##使用monocle选择的高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
p3 <- plot_ordering_genes(mycds)
p3
#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)
#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
plot1
ggsave("Epi cluster/serrated_epi pseudotime/State.pdf", plot = plot1, width = 6, height = 5)
ggsave("Epi cluster/serrated_epi pseudotime/State.png", plot = plot1, width = 6, height = 5)
##Cluster轨迹分布图
plot2 <- plot_cell_trajectory(mycds, color_by = "celltypes4")
plot2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
plot2 <-plot_cell_trajectory(mycds, color_by = "celltypes4") + 
  scale_color_manual(breaks = c("Transit amplifying cells", 'Intermediate',"Epi-1", "Epi-2"), 
                     values=c("#CAB2D6", 
                              "#FDBF6F","#33A02C","#E31A1C")) +
  theme(legend.position = "right", axis.text = element_text(size = 16), axis.title = element_text(size = 16))
plot2 
ggsave("Epi cluster/Epi2 pseudotime/Cluster20211203 V2.pdf", plot = plot2, width = 4, height = 3)
ggsave("Epi cluster/serrated_epi pseudotime/TAINTEPI12 V2 20211206.pdf", plot = plot2, width = 5*1.5, height = 4.5*1.5)
##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")+
  theme(legend.position = "right", axis.text = element_text(size = 16), axis.title = element_text(size = 16))

plot3
ggsave("Epi cluster/serrated_epi pseudotime/Pseudotime20210726 V2 20211206.pdf", plot = plot3, width =5*1.5/1.15, height = 4.5*1.5)
ggsave("Epi cluster/serrated_epi pseudotime/Pseudotime.png", plot = plot3, width = 6/1.5, height = 4.5/1.5)


plot_cell_trajectory(mycds, color_by = "celltypes4") +facet_wrap(~celltypes4, ncol = 3)
plot_cell_trajectory(mycds, color_by = "State") + facet_wrap(~State, ncol = 3)
mycds <- orderCells(mycds, root_state = 3)




## Figure 3D ## 
library(Seurat)
load("i:/genomicdata/jiaoda/raw/seurat_epi.Rdata")
load("i:/genomicdata/external/ZYJ/12/Epi_metadata20211125V2.Rdat")
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

seurat@meta.data <- Epi_metadata
mtx <- AverageExpression(object = seurat, group.by = "celltypes4")[[1]]
# plot.gene <- c(paste0("ANXA", seq(1:10)), "S100P", "S100A2", "S100A4", "S100A6", "S100A8", "S100A11")
plot.gene <- c("ANXA10", "ASXL1", "CFTR", "DOT1L", "HIC1", "INO80", "KLF3", "MCM3AP", "MCM8", "PDLIM2", "POLD1", "TP53BP1", "WNK2", "WRN")
plot.celltype <- c("Epi-1", "Epi-2", "Transit amplifying cells", "Goblet cells", "Intermediate", "BEST4+ colonocytes", "Mature colonocytes")
plot.data <- mtx[plot.gene, plot.celltype]
plot.data <- standarize.fun(indata = plot.data, halfwidth = 2)
annCol <- data.frame(celltypes = factor(x = colnames(plot.data), levels = colnames(plot.data)),
                     row.names = colnames(plot.data),
                     stringsAsFactors = F)
annColors <- list()
annColors[["celltypes"]] = c("Epi-2" = "#E31A1C",
                             # "Epi-3" = "#E10C0D",
                             "Epi-1" = "#33A02C",
                             "Transit amplifying cells" = "#C9B1D5",
                             "Goblet cells" = "#FB9998",
                             "Intermediate" = "#FDBE6E",
                             "BEST4+ colonocytes" = "#A6CEE3",
                             "Mature colonocytes" = "#B2DF8A",
                             "Tuft cells" = "#FF7F00",
                             "Enteroendocrine cells" = "#6A3D9A") 

pdf("i:/genomicdata/jiaoda/plot/Heatmap/Epi1203-1.pdf", width = 587/100, height = 412/100)
pheatmap(plot.data,
         cluster_cols = F, cluster_rows = F,
         border_color = "white",
         color = colorRampPalette(colors = c("#6353A3","white","#9A0038"))(100),
         annotation_col = annCol,
         annotation_colors = annColors,
         cellwidth = 12, cellheight = 12,
         # angle_col = "45",
         treeheight_row = 20)
invisible(dev.off())



## Figure 3E ## 
x <- read.csv("Epi cell_identify/EpiSL cell_identify/Epi1 vs 2.diff_genes_wilcox20211125V2.csv", row.names = 1,check.names = F,stringsAsFactors = F,header = T)
x$label<- rownames(x)
head(x)
colnames(x)[2] <- "logFC"
colnames(x)[5] <- "padj"

# 突出展示感兴趣的基因
selectedGeneID <- c('TFF1',
                    'LYZ',
                    'CEACAM6',
                    'S100A11',
                    'CEACAM5',
                     'S100A6',
                    'TSPAN8',
                    'S100P',
                    'DUOX2',
                    'SLC40A1',
                    'DUOXA2',
                    'ANXA1',
                    'TFF2',
                    'JUN',
                    'SRSF7',
                    'MKI67',
                    'SNHG12',
                    'HES1',
                    'FOSB',
                    'FABP5',
                    'SOX4',
                    'EGR1',
                    'HMGB1',
                    'ELF3',
                    'MUC2',
                    'ANXA10',
                    'NQO1',
                    'GSTP1',
                    'MUC5B',
                    'GPX4',
                    'ATP5ME',
                    'CDK1',
                    'JUNB',
                    'ATP5MC1',
                    'MUC4',
                    'DLL4',
                    'MYCL',
                    'MUC5AC',
                    'BRCA1',
                    'DUSP2',
                    'SOD2',
                    'SRSF1',
                    'OLFM4'
)

selectgenes  <- x[selectedGeneID,]

plot_mode <- "classic" #经典

logFCcut <- 0.5 #log2-foldchange
pvalCut <- 0.05 #P.value
adjPcut <- 0.05 #adj.P.value


#置x，y軸的最大最小
xmin <- (range(x$logFC)[1]- (range(x$logFC)[1]+ 10))
xmax <- (range(x$logFC)[1]+ (10-range(x$logFC)[1]))
ymin <- 0
ymax <- 250

# 簡單的setting for color
x$color_transparent <- ifelse((x$padj < adjPcut & x$logFC > logFCcut), heatmap.L.BlYlRd[5], ifelse((x$padj < adjPcut & x$logFC < -logFCcut), heatmap.L.BlYlRd[1],heatmap.L.BlYlRd[3]))
# 簡單的setting for size
size <- ifelse((x$padj < adjPcut & abs(x$logFC) > logFCcut), 4, 2)

p1 <- ggplot(data=x, aes(logFC, -log10(padj), label = label)) +
  geom_point(alpha = 0.8, size = size, colour = x$color_transparent) +
  
  labs(x="log2FoldChange", y="log10FDR", title="") + 
  ylim(c(ymin,ymax)) + 
  scale_x_continuous(
    breaks = c(-2.5,-2,-1.5, -1, -logFCcut, 0, logFCcut, 1, 1.5,2), #刻度线的位置
    labels = c(-2.5,-2,-1.5, -1, -logFCcut, 0, logFCcut, 1, 1.5,2),
    limits = c(-2.75, 2) #x轴范围
  ) +
  scale_y_continuous(
    breaks = c(0,50,100,150,200,250),
    labels = c(0,50,100,150,200,250),
    limits = c(0,270)
  ) +
  #或用下面这行???
  #xlim(c(xmin, xmax)) + 
  
  #画阈值分界线
  geom_vline(xintercept = c(-logFCcut, logFCcut), color="grey40", 
             linetype="longdash", lwd = 0.5) + #虚线的形状和粗细
  geom_hline(yintercept = -log10(pvalCut), color="grey40", 
             linetype="longdash", lwd = 0.5) +
  
  theme_bw(base_size = 12#, base_family = "Times" #修改字体
  ) +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        panel.background = element_blank())
p1

p2 <- p1 + 
  # 在感兴趣的基因外面画个黑色圈
  geom_point(data = selectgenes, alpha = 1, size = 5, shape = 1, #size白圈
             stroke = 1, #圈粗???
             color = heatmap.L.BlYlRd[5]) +
  
  # 显示感兴趣的基因的基因名
  #scale_color_manual(values = mycol) + 
  geom_text_repel(data = selectgenes, 
                  show.legend = FALSE, 
                  size = 4, box.padding = unit(0.35, "lines"), 
                  point.padding = unit(0.3, "lines")) +
  guides(color=guide_legend(title = NULL)) 

p2 <- p2 + coord_flip()
p2
ggsave("volcano 20211205V2.pdf", width = 6,height = 6)




## Figure 3F ## 
# load R package
library(GSVA)
library(survival)
library(survminer)
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)
library(ggplot2)
library(ggpubr)
library(gplots)
library(pheatmap)
library(ComplexHeatmap)
library(viridis)
library(clusterProfiler)

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

# loac color
jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")
load("cell.col.RData")
cellLevel <- c('BEST4+ colonocytes','Mature colonocytes','Goblet cells','Intermediate',
               'Transit amplifying cells','Tuft cells',
               'Enteroendocrine cells',
               'Epi-1','Epi-2')

# load data
load("C2gsva_epi.Rdat")
epi.c2 <- dat2
load("Hallmarkgsva_epi.Rdat")
epi.h <- dat2
rm(dat2); gc()
epires <- readRDS("Epires.rds")

# draw gsva heatmap
hmat <- epi.h[,c("celltypes",
                 "HALLMARK_MYC_TARGETS_V1",
                 "HALLMARK_MYC_TARGETS_V2",
                 "HALLMARK_G2M_CHECKPOINT",
                 "HALLMARK_NOTCH_SIGNALING",
                 "HALLMARK_E2F_TARGETS")]
rownames(hmat) <- epi.h$X
hmat$celltypes <- factor(hmat$celltypes,levels = cellLevel)
hmat <- hmat[order(hmat$celltypes),]

c2mat <- epi.c2[,c("celltypes",
                   "BIOCARTA_DNAFRAGMENT_PATHWAY",
                   "LEE_METASTASIS_AND_RNA_PROCESSING_UP",
                   "GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_TAN_DN",
                   "COLLER_MYC_TARGETS_UP",
                   "WONG_EMBRYONIC_STEM_CELL_CORE",
                   "DANG_MYC_TARGETS_UP")]
rownames(c2mat) <- epi.c2$X

c2mat$celltypes <- factor(c2mat$celltypes,levels = cellLevel)
c2mat <- c2mat[order(c2mat$celltypes),]

kmat <- epi.c2[,c("celltypes",
                  "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
                  "KEGG_ARGININE_AND_PROLINE_METABOLISM",
                  "KEGG_FATTY_ACID_METABOLISM",
                  "KEGG_CITRATE_CYCLE_TCA_CYCLE",
                  "KEGG_BUTANOATE_METABOLISM",
                  "REACTOME_BETA_OXIDATION_OF_LAUROYL_COA_TO_DECANOYL_COA_COA",
                  "REACTOME_BETA_OXIDATION_OF_HEXANOYL_COA_TO_BUTANOYL_COA")]
rownames(kmat) <- epi.c2$X
kmat$celltypes <- factor(kmat$celltypes,levels = cellLevel)
kmat <- kmat[order(kmat$celltypes),]

annCol <- hmat[,1,drop = F]
annColors <- list(celltypes = cell.col)

annCol2 <- annCol[-which(annCol$celltypes %in% c("Tuft cells","Enteroendocrine cells","Epi-1")),,drop = F]
annColors2 <- list(celltypes = annColors$celltypes[c(1,2,3,4,5,6,8)])
indata1 <- hmat[,-1]
indata1 <- standarize.fun(t(indata1), halfwidth = 1)
rownames(indata1) <- gsub("_"," ",rownames(indata1))
rownames(indata1) <- tolower(rownames(indata1))
hm1 <- pheatmap(as.matrix(indata1[,rownames(annCol2)]),
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                border_color = NA,
                color = greenred(64),
                annotation_col = annCol2,
                annotation_colors = annColors2,
                gaps_col = c(2634,3654,5311),
                cellwidth = 0.05,
                cellheight = 12,
                name = "HALLMARK")

indata2 <- c2mat[,-1]
indata2 <- standarize.fun(t(indata2), halfwidth = 1)
rownames(indata2) <- gsub("_"," ",rownames(indata2))
rownames(indata2) <- tolower(rownames(indata2))
hm2 <- pheatmap(as.matrix(indata2[,rownames(annCol2)]),
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                border_color = NA,
                color = inferno(64),
                gaps_col = c(2634,3654,5311),
                cellwidth = 0.05,
                cellheight = 12,
                name = "C2")

indata3 <- kmat[,-1]
tmp <- rbind.data.frame(epires$KEGG[c("Fatty acid biosynthesis",
                                      "One carbon pool by folate"),],
                        epires$REACTOME[c("Glucose metabolism",
                                          "Metabolism of folate and pterines",
                                          "Foxo mediated transcription of oxidative stress metabolic and neuronal genes"),])
colnames(tmp) <- gsub(".","-",colnames(tmp),fixed = T)
indata3 <- cbind.data.frame(indata3,t(tmp[,rownames(indata3)]))
colnames(indata3) <- gsub("_"," ",colnames(indata3))
colnames(indata3) <- tolower(colnames(indata3))
indata3 <- standarize.fun(t(indata3), halfwidth = 1)
hm3 <- pheatmap(as.matrix(indata3[,rownames(annCol2)]),
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                border_color = NA,
                color = viridis(64),
                gaps_col = c(2634,3654,5311),
                cellwidth = 0.05,
                cellheight = 12,
                name = "KEGG")

pdf("heatmap of three pathway categories2.pdf", width = 14,height = 8)
draw(hm1 %v% hm2 %v% hm3, annotation_legend_side = "bottom",heatmap_legend_side = "bottom", use_raster = F)
invisible(dev.off())



## Figure 3G ## 
# boxplot for emtab7960
emtab7960.expr <- read.table("alldata.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
emtab7960.expr <- emtab7960.expr[,-2]
emtab7960.expr <- apply(emtab7960.expr[,setdiff(colnames(emtab7960.expr), "Gene")], 2, function(x) tapply(x, INDEX=factor(emtab7960.expr$Gene), FUN=mean, na.rm=TRUE))
emtab7960.expr <- as.data.frame(emtab7960.expr)
emtab7960.surv <- read.table("E7960_clincal.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
emtab7960.expr <- emtab7960.expr[,rownames(emtab7960.surv)]
emtab7960.surv$group2 <- ifelse(emtab7960.surv$group == "ND","ND","HD_LD")
epi1.degs <- read.csv("Epi1 vs normal.diff_genes_wilcox.csv",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
epi2.degs <- read.csv("Epi2 vs normal.diff_genes_wilcox.csv",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
epi3.degs <- read.csv("Epi3 vs normal.diff_genes_wilcox.csv",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
DEfiles <- c("Epi1 vs normal.diff_genes_wilcox.csv",
             "Epi2 vs normal.diff_genes_wilcox.csv",
             "Epi3 vs normal.diff_genes_wilcox.csv")
dirct <- "up"
p.cutoff <- 0.05
p.adj.cutoff <- 0.05
n.marker <- 50
genelist <- c()
for (filek in DEfiles) {
  DEres <- read.csv( filek, header=TRUE, row.names=1, stringsAsFactors=FALSE,check.names = F)
  #rownames(DEres) <- toupper(rownames(DEres))
  if (dirct == "up") {
    genelist <- c( genelist, rownames(DEres[!is.na(DEres$p_val_adj) & DEres$p_val < p.cutoff & DEres$p_val_adj < p.adj.cutoff & !is.na(DEres$avg_logFC) & DEres$avg_logFC > 0, ]) )
  }
  if (dirct == "down") {
    genelist <- c( genelist, rownames(DEres[!is.na(DEres$p_val_adj) & DEres$p_val < p.cutoff & DEres$p_val_adj < p.adj.cutoff & !is.na(DEres$avg_logFC) & DEres$avg_logFC < 0, ]) )
  }
}
unqlist <- setdiff(genelist,genelist[duplicated(genelist)])

marker <- list()
for (filek in DEfiles) {
  DEres <- read.csv( filek, header=TRUE, row.names=1, stringsAsFactors=FALSE,check.names = F)
  #rownames(DEres) <- toupper(rownames(DEres))
  if(dirct == "up") {
    outk <- intersect( unqlist, rownames(DEres[!is.na(DEres$p_val_adj) & DEres$p_val < p.cutoff & DEres$p_val_adj < p.adj.cutoff & !is.na(DEres$avg_logFC) & DEres$avg_logFC > 0, ]) )
    outk <- DEres[outk,]
    outk <- outk[order(outk$avg_logFC, decreasing = TRUE),]
    
    if(nrow(outk) > n.marker) {
      marker[[filek]] <- outk[1:n.marker,]
    } else {
      marker[[filek]] <- outk
    }
    marker$dirct <- "up"
  }
  if(dirct == "down") {
    outk <- intersect( unqlist, rownames(DEres[!is.na(DEres$p_val_adj) & DEres$p_val < p.cutoff & DEres$p_val_adj < p.adj.cutoff & !is.na(DEres$avg_logFC) & DEres$avg_logFC < 0, ]) )
    outk <- DEres[outk,]
    outk <- outk[order(outk$avg_logFC, decreasing = FALSE),]
    
    if(nrow(outk) > n.marker) {
      marker[[filek]] <- outk[1:n.marker,]
    } else {
      marker[[filek]] <- outk
    }
    marker$dirct <- "down"
  }
}

# generate templates for calculate enrichment score
templates <- NULL
for (filek in DEfiles) {
  tmp <- data.frame(probe = rownames(marker[[filek]]),
                    class = sub(" vs normal.diff_genes_wilcox.csv","",filek),
                    dirct = marker$dirct,
                    stringsAsFactors = FALSE)
  templates <- rbind.data.frame(templates, tmp, stringsAsFactors = FALSE)
}
#write.table(templates, file=paste0(dirct,"regulated_marker_templates_epi.txt"), row.names = FALSE, sep = "\t", quote = FALSE)

# generate enrichment score
signature <- list()
signature[["Epi1"]] <- templates[which(templates$class == "Epi1"),"probe"]
signature[["Epi2"]] <- templates[which(templates$class == "Epi2"),"probe"]
signature[["Epi3"]] <- templates[which(templates$class == "Epi3"),"probe"]
epi.score <- gsva(as.matrix(log2(emtab7960.expr + 1)),
                  signature,
                  method = "ssgsea")
emtab7960.surv$Epi1 <- as.numeric(epi.score[1,rownames(emtab7960.surv)])
emtab7960.surv$Epi2 <- as.numeric(epi.score[2,rownames(emtab7960.surv)])
emtab7960.surv$Epi3 <- as.numeric(epi.score[3,rownames(emtab7960.surv)])
boxplot(Epi1~group2,emtab7960.surv)
boxplot(Epi2~group2,emtab7960.surv)
boxplot(Epi3~group2,emtab7960.surv)

wilcox.test(Epi1~group2,emtab7960.surv)
wilcox.test(Epi2~group2,emtab7960.surv)
wilcox.test(Epi3~group2,emtab7960.surv)
emtab7960.surv$group2 <- factor(emtab7960.surv$group2,levels = c("ND","HD_LD"))
ggplot(data = emtab7960.surv,aes(x = group2, #分组列名
                                 y = Epi2, #连续变量列名
                                 fill = group2))+ #按分组填充颜???
  scale_fill_manual(values = jco[1:2]) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="white") +
  
  geom_point(shape = 21, size=2, # 点的性状和大???
             position = position_jitterdodge(), # 让点散开
             aes(color=group2), alpha = 1) +
  scale_color_manual(values =jco[1:2]) + #用自定义颜色填充
  geom_boxplot(notch = FALSE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  scale_x_discrete(labels= c("SSL\n(n=22)","SSL\nwith dysplasia\n(n=20)")) + 
  theme_bw() + 
  ylab("Epi-2 signature\nssgsea [z-scored") +
  xlab("") +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  stat_compare_means(method = "wilcox",label.y = min(emtab7960.surv$Epi2))
ggsave("boxplot for epi2.pdf",width = 2.6,height = 4)

ggplot(data = emtab7960.surv,aes(x = group2, #分组列名
                                 y = Epi3, #连续变量列名
                                 fill = group2))+ #按分组填充颜???
  scale_fill_manual(values = jco[1:2]) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="white") +
  
  geom_point(shape = 21, size=2, # 点的性状和大???
             position = position_jitterdodge(), # 让点散开
             aes(color=group2), alpha = 1) +
  scale_color_manual(values =jco[1:2]) + #用自定义颜色填充
  geom_boxplot(notch = FALSE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 0.7) +
  scale_x_discrete(labels= c("SSL\n(n=22)","SSL\nwith dysplasia\n(n=20)")) + 
  theme_bw() + 
  ylab("Epi-3 signature\nssgsea [z-scored") +
  xlab("") +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  stat_compare_means(method = "wilcox",label.y = min(emtab7960.surv$Epi3))
ggsave("boxplot for epi3.pdf",width = 3,height = 4)

# 
msigdb <- read.gmt("msigdb.v7.4.symbols.gmt")
#msigdb <- msigdb[which(msigdb$term %in% c(rownames(indata1),rownames(indata2),rownames(indata3))),]
msigdb <- msigdb[which(msigdb$term %in% c(rownames(indata1),rownames(indata2)[-c(1,7)])),]
msigdb$term <- as.character(msigdb$term)
msigdb$gene <- as.character(msigdb$gene)
cell.type <- unique(msigdb$term)
pathway <- list()
for (i in cell.type) {
  pathway[[i]] <- msigdb[which(msigdb$term == i),"gene"]
}

emtab7960.ssgsea <- gsva(as.matrix(log2(emtab7960.expr + 1)),
                         pathway,
                         method = "ssgsea")
pval <- c()
for (i in rownames(emtab7960.ssgsea)) {
  tmp <- data.frame(score = as.numeric(emtab7960.ssgsea[i,]),
                    group = emtab7960.surv$group2)
  wt <- wilcox.test(tmp$score~tmp$group)
  pval <- c(pval,wt$p.value)
}
pval <- format(round(pval,3),nsmall = 2)
pval <- ifelse(pval == "0.000","P<0.001",paste0("P=",pval))
names(pval) <- rownames(emtab7960.ssgsea)

plotdata <- standarize.fun(emtab7960.ssgsea,halfwidth = 1)
plotdata1 <- plotdata[rownames(indata1),]; rownames(plotdata1) <- paste0(rownames(plotdata1)," ",pval[rownames(plotdata1)])
plotdata2 <- plotdata[rownames(indata2)[-c(1,7)],]; rownames(plotdata2) <- paste0(rownames(plotdata2)," ",pval[rownames(plotdata2)])
#plotdata3 <- plotdata[rownames(indata3),]; rownames(plotdata3) <- paste0(rownames(plotdata3)," ",pval[rownames(plotdata3)])

annCol.emtab <- emtab7960.surv[,"group2",drop = F]
annCol.emtab$group2 <- ifelse(annCol.emtab$group2 == "ND","SSL","SSL with dysplasia")
annCol.emtab <- annCol.emtab[order(annCol.emtab$group2),,drop = F]
annColors.emtab <- list("group2" = c("SSL" = jco[1],"SSL with dysplasia" = jco[2]))

hm1 <- pheatmap(as.matrix(plotdata1[,rownames(annCol.emtab)]),
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                border_color = NA,
                color = greenred(64),
                annotation_col = annCol.emtab,
                annotation_colors = annColors.emtab,
                cellwidth = 8,
                cellheight = 12,
                name = "HALLMARK")

hm2 <- pheatmap(as.matrix(plotdata2[,rownames(annCol.emtab)]),
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                border_color = NA,
                color = inferno(64),
                cellwidth = 8,
                cellheight = 12,
                name = "C2")

# hm3 <- pheatmap(as.matrix(plotdata3[,rownames(annCol.emtab)]),
#                 cluster_rows = F,
#                 cluster_cols = F,
#                 show_rownames = T,
#                 show_colnames = F,
#                 border_color = NA,
#                 color = viridis(64),
#                 cellwidth = 8,
#                 cellheight = 12,
#                 name = "KEGG")

pdf("heatmap of three pathway categories in EMTAB7960.pdf", width = 14,height = 6)
draw(hm1 %v% hm2, annotation_legend_side = "bottom",heatmap_legend_side = "bottom")
invisible(dev.off())



## Figure 3H ## 
library(Seurat)
load("i:/genomicdata/jiaoda/raw/seurat_epi.Rdata")
load("i:/genomicdata/external/ZYJ/12/Epi_metadata20211125V2.Rdat")
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

seurat@meta.data <- Epi_metadata
mtx <- AverageExpression(object = seurat, group.by = "celltypes4")[[1]]
plot.gene <- c("JAG1", "DLL1", "DLL4", "NOTCH1", "MYC", "MYCL", "HES1")
plot.celltype <- c("Epi-1", "Epi-2", "Transit amplifying cells", "Goblet cells", "Intermediate", "BEST4+ colonocytes", "Mature colonocytes")
plot.data <- mtx[plot.gene, plot.celltype]
plot.data <- standarize.fun(indata = plot.data, halfwidth = 2)
annCol <- data.frame(celltypes = factor(x = colnames(plot.data), levels = colnames(plot.data)),
                     row.names = colnames(plot.data),
                     stringsAsFactors = F)
annColors <- list()
annColors[["celltypes"]] = c("Epi-2" = "#E31A1C",
                             # "Epi-3" = "#E10C0D",
                             "Epi-1" = "#33A02C",
                             "Transit amplifying cells" = "#C9B1D5",
                             "Goblet cells" = "#FB9998",
                             "Intermediate" = "#FDBE6E",
                             "BEST4+ colonocytes" = "#A6CEE3",
                             "Mature colonocytes" = "#B2DF8A",
                             "Tuft cells" = "#FF7F00",
                             "Enteroendocrine cells" = "#6A3D9A") 

pdf("i:/genomicdata/jiaoda/plot/Heatmap/Epi1203-1.pdf", width = 587/100, height = 412/100)
pheatmap(plot.data,
         cluster_cols = F, cluster_rows = F,
         border_color = "white",
         color = colorRampPalette(c("#75A3D0", "white", "#FF0303"))(64),
         annotation_col = annCol,
         annotation_colors = annColors,
         cellwidth = 12, cellheight = 12,
         # angle_col = "45",
         treeheight_row = 20)
invisible(dev.off())



## Figure 3J ##
plot1<-VlnPlot(seurat, features = c('LYZ'), pt.size=0, group.by='celltypes4', ncol=3)+geom_boxplot(width=0.2,col="black",fill="#E1E3E6", outlier.size = 0.5)
plot1
ggsave("备用/Figure S2/Vlnplot LYZ Epi2.pdf", plot = plot1, width = 15, height = 4)


