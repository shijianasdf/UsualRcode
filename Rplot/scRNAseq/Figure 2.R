##### Figure plotting #####
## Figure 2A left ## 
library(Seurat)
library(ggplot2)
library(dplyr)
###读入
setwd('D:/2021scRNA-seq')
load("seurat_epi.Rdat")
###统计各cluster比例
prop.group<-as.data.frame(table(seurat$group))
table(seurat$new.patient)
prop.cluster<-as.data.frame(table(seurat@meta.data$seurat_clusters,seurat$group))
colnames(prop.cluster)<-c('cluster','group','freq.cluster')
prop_data <- merge(prop.group,prop.cluster,by.x=c("Var1"),by.y=c("group"))
prop_data$clusterprop<-100*prop_data$freq.cluster/prop_data$Freq
arrange(prop_data,cluster)

seurat  <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 4000)
top10<-head(VariableFeatures(seurat),10)
top10
scale.genes <-  rownames(seurat)
seurat <- ScaleData(seurat, features = scale.genes)
seurat <- RunPCA(seurat , features = VariableFeatures(object = seurat ))
plot2 <- ElbowPlot(seurat, ndims=20, reduction="pca") 
plot2
pc.num=1:11
seurat <- FindNeighbors(seurat, dims = pc.num)
seurat <- FindClusters(seurat, resolution = 0.5)
table(seurat@meta.data$seurat_clusters)
metadata <- seurat@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'Epi cluster/cell_cluster nfeatures4000 20210325.csv',row.names = F)

seurat <- RunUMAP(seurat, n.neighbors = 50,n.epochs=50,dims = pc.num)
embed_umap <- Embeddings(seurat, 'umap')
write.csv(embed_umap,'Epi cluster/embed_umap n.neighbors = 50,n.epochs=50 20210325.csv') 
plot2 = DimPlot(seurat,label = T, reduction = "umap") 
plot2
ggsave("Epi cluster/UMAP 20210411.pdf", plot = plot2, width = 8, height = 7)
ggsave("Epi cluster/UMAP n.neighbors = 50,n.epochs=50 20210323.png", plot = plot2, width = 8, height = 7)
diff.wilcox = FindAllMarkers(seurat)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(all.markers, "Epi cell_identify/diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "Epi cell_identify/top10_diff_genes_wilcox.csv", row.names = F)

Idents(seurat) <- "Cell_subtype" 
DimPlot(seurat, reduction = "umap",group.by = "group",label = T)+ theme_bw()
DimPlot(seurat, reduction = "umap",group.by = "Cell_subtype",label = T)+ theme_bw()
DimPlot(seurat, reduction = "umap",group.by = "patient",label = F)

#根据以上注释结果对细胞群添加注释信息
seurat@meta.data$celltypes2<-seurat@meta.data$seurat_clusters
Idents(seurat) <- "celltypes2" 
seurat<- RenameIdents(seurat, '0'='Epi-SL','2'='Epi-SL','12'='Epi-SL','6'='Epi-SL',
                      '7'='Epi-SL','10'='Epi-SL','13'='Epi-SL','3'='Transit amplifying cells',
                      '4'='Goblet cells','5'='Goblet cells','11'='Goblet cells','8'='Intermediate',
                      '9'='BEST4+ colonocytes','1'='Mature colonocytes','14'='Tuft cells',
                      '15'='Enteroendocrine cells')
seurat$celltypes2<-Idents(seurat)
plot2<-DimPlot(seurat, reduction = "umap",cols=c("#074083","#CAB2D6", "#FB9A99",
                                                 "#FDBF6F", "#A6CEE3", "#B2DF8A" ,"#FF7F00" ,"#6A3D9A"),
               group.by ='celltypes2',label = T)
plot2
ggsave("Epi cluster/UMAP 20210726.pdf", plot = plot2, width = 7.5, height = 5)



## Figure 2A right ## 

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
Epi_metadata<-seurat@meta.data
dat <- Epi_metadata[,c("new.patient","group","celltypes2")]
dat$new.patient <- factor(dat$new.patient,levels = rev(c("NC1","NC2","NC3","HP1","HP2","HP3","HP4","SSL1","SSL2","SSL3","SSL4","SSL5","SSLD","TSA1","TSA2","TSA3","TSA4","TSA5")))
dat$celltypes <- factor(dat$celltypes2,levels = rev(c('BEST4+ colonocytes','Mature colonocytes','Goblet cells','Intermediate',
                                                      'Transit amplifying cells','Epi-SL','Tuft cells',
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
              "Epi-SL" = "#074083",
              "Tuft cells" = mycol[9],
              "Enteroendocrine cells" = mycol[10])
save(cell.col,file = "cell.col.RData")

dat.bar <- as.data.frame(table(dat$new.patient,dat$celltypes))
dat.bar <- dat.bar%>% group_by(Var1) %>% 
  mutate(Pct = Freq / sum(Freq) * 100) %>% as.data.frame()
colnames(dat.bar) <- c("patient","celltype","count","pct")

# left barplot

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
ggsave(filename = "distribution of epi cells2 20211125V2.pdf", width = 5,height = 3)
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
ggsave(filename = "distribution of epi cells2.pdf", width = 5,height = 3)



## Figure 2B ## 
load("seurat_epi.Rdat")
marker <- list("Mature colonocytes" = c("GUCA2B","SLC26A3","CA1"),
               "BEST4+ colonocytes" = c("BEST4","OTOP2","LYPD8","CA4","CA7","MT1H"),
               "Goblet cells" = c("TFF3","SPINK1","SPINK4","REG4","AGR2","MUC2","WFDC2","CLCA1"),
               "Transit amplifying cells" = c("MKI67","PCNA","LGR5","ASCL2","OLFM4","CDX2","ATP5MC1","ATP5MC3","LYZ"),
               "Tuft cells" = c("HPGDS","PTGS1"),
               "Enteroendocrine cells" = c("SCGN","PCSK1N"))
module.gene <- sapply(marker,length)
cellLevel <- c('BEST4+ colonocytes','Mature colonocytes','Goblet cells','Intermediate',
               'Transit amplifying cells','Tuft cells',
               'Enteroendocrine cells',
               'Epi-SL')
marker <- rev(unlist(marker[cellLevel]))
annotation <- data.frame(cell = seurat$orig.ident,
                         cluster = as.numeric(as.character(seurat$seurat_clusters)),
                         celltype = seurat$celltypes2)

# extract averaged expression
marker.expr.celltype <- AverageExpression(seurat,assays = "RNA",features = marker,group.by = "celltypes2",verbose = TRUE)
marker.expr.celltype <- as.data.frame(marker.expr.celltype$RNA)
marker.expr.celltype <- marker.expr.celltype[,cellLevel]
marker.expr.celltype$gene <- rownames(marker.expr.celltype)


# calculate expression fraction
fraction.celltype <- matrix(0,nrow = length(marker),ncol = length(cellLevel),byrow = T,
                            dimnames = list(marker,cellLevel))

for (i in colnames(fraction.celltype)) {
  message("--",i)
  tmp <- FetchData(object = subset(seurat, celltypes2 == i), vars = marker)
  for (j in marker) {
    fraction.celltype[j,i] <- sum(as.numeric(tmp[,j]) > 0) /nrow(tmp)
  }
}
fraction.celltype <- as.data.frame(fraction.celltype); fraction.celltype$gene = rownames(fraction.celltype)
fraction.celltype <- fraction.celltype[,colnames(marker.expr.celltype)]
identical(colnames(fraction.celltype),colnames(marker.expr.celltype))

# generate cell type dotplot
fraction.celltype.long <- gather(fraction.celltype, celltype, pct, `BEST4+ colonocytes`:`Epi-SL`, factor_key=TRUE)
expr.celltype.long <- gather(marker.expr.celltype, celltype, expr, `BEST4+ colonocytes`:`Epi-SL`, factor_key=TRUE)

tmp <- cbind.data.frame(fraction.celltype.long,expr = log2(expr.celltype.long$expr + 1))
tmp$celltype <- factor(tmp$celltype,levels = cellLevel)
tmp$gene <- factor(tmp$gene,levels = marker)
ggplot(tmp,aes(x = celltype, y = gene,size = pct)) +
  geom_point(shape = 21, aes(fill = expr), position = position_dodge(0), stroke = 1) +
  xlab(NULL) + ylab(NULL) +
  labs(size = "Proportion of\nexpressing\ncells", fill = "Normalized\ngene exp.") +
  scale_x_discrete(position = "top") + 
  scale_size_continuous(range = c(1,6)) + 
  #scale_fill_gradient2(low = inferno(3)[1],mid = inferno(3)[2],high = inferno(3)[3],midpoint = sum(range(log2(expr.celltype.long$expr + 1)))/2) +
  #scale_fill_gradient2(low = brewer.pal(11,"Spectral")[11],mid = "white",high = brewer.pal(11,"Spectral")[1],midpoint = sum(range(log2(expr.celltype.long$expr + 1)))/2) +
  scale_fill_gradient2(low = brewer.pal(11,"Spectral")[11],mid = brewer.pal(11,"Spectral")[6],high = brewer.pal(11,"Spectral")[1],midpoint = sum(range(log2(expr.celltype.long$expr + 1)))/2) +
  theme_bw() + 
  theme(legend.position = "right",
        legend.box = "vertical",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 10, colour = "black",angle = 90, vjust = 0.5, hjust=0),
        axis.text.y = element_text(size = 10, colour = "black"),
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0.1,"in")) +
  geom_hline(yintercept = cumsum(rev(na.omit(module.gene[cellLevel])))[-length(module.gene)] + 0.5, color="grey40", linetype="longdash", size=0.5) +
  geom_vline(xintercept = 7.5, color="grey40", linetype="longdash", size=0.5)
ggsave(filename = "dotplot for avgExpr and pct among cell types20211125V2.pdf", width = 5,height = 8)



## Figure 2C ## 
x <- read.csv("Epi cell_identify/EpiSL(V2) vs normal coloncyte.diff_genes_wilcox.csv", row.names = 1,check.names = F,stringsAsFactors = F,header = T)
x$label<- rownames(x)
head(x)
colnames(x)[2] <- "logFC"
colnames(x)[5] <- "padj"

# 突出展示感兴趣的基因
selectedGeneID <- c('PKM',
                    'LYZ',
                    'GSTP1',
                    'GOLM1',
                    'HES1',
                    'OLFM4',
                    'S100A11',
                    'ANXA2',
                    'GPX2',
                    'ENO1',
                    'JAG1',
                    'SOX4',
                    'ATP5ME',
                    'LCN2',
                    'MKI67',
                    'FOXQ1',
                    'NQO1',
                    'S100P',
                    'DUOX2',
                    'MYC',
                    'DUOXA2',
                    'SOX9','SERPINB6'
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
    breaks = c(-1.5, -1, -logFCcut, 0, logFCcut, 1, 1.5,2), #刻度线的位置
    labels = c(-1.5, -1, -logFCcut, 0, logFCcut, 1, 1.5,2),
    limits = c(-1.5, 2.1) #x轴范围
  ) +
  scale_y_continuous(
    breaks = c(0,50,100,150),
    labels = c(0,50,100,150),
    limits = c(0,180)
  ) +

  
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
ggsave("volcano 20211203V2.pdf", width = 6,height = 5.5)



## Figure 2D ## 
#####与正常细胞比做GO富集
EpiSL.de.markers2 <- FindMarkers(seurat, ident.1 = "Epi-SL", ident.2 = c('Mature colonocytes','BEST4+ colonocytes','Goblet cells'))
write.csv(EpiSL.de.markers2, "Epi cell_identify/EpiSL(V2) vs normal coloncyte.diff_genes_wilcox.csv", row.names = T)

EpiSLG <- read.csv(file = "Epi cell_identify/GO/EpiSL2全.csv", header = T)
genes <- EpiSLG[[1]]

#GO富集分析
library(org.Hs.eg.db)
library(clusterProfiler)
#对于加载的注释库的使用，以上述为例，就直接在 OrgDb 中指定人（org.Hs.eg.db）或绵羊（sheep）
enrich.go <- enrichGO(gene = genes,  #基因列表文件中的基因名称
                      OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库，示例物种是绵羊（sheep）
                      keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                      ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'fdr',  #指定 p 值校正方法
                      pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                      qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                      readable = FALSE)

#例如上述指定 ALL 同时计算 BP、MF、CC，这里将它们作个拆分后输出
BP <- enrich.go[enrich.go$ONTOLOGY=='BP', ]
CC <- enrich.go[enrich.go$ONTOLOGY=='CC', ]
MF <- enrich.go[enrich.go$ONTOLOGY=='MF', ]

write.table(as.data.frame(BP), 'Epi cell_identify\\EPISL(V2)DEG\\go.BP.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(as.data.frame(CC), 'Epi cell_identify\\EPISL(V2)DEG\\go.CC.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(as.data.frame(MF), 'Epi cell_identify\\EPISL(V2)DEG\\go.MF.txt', sep = '\t', row.names = FALSE, quote = FALSE)

gene = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gene)

ego <- enrichKEGG(
  gene          = gene$ENTREZID,
  keyType     = "kegg",
  organism   = 'hsa',
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.1
)
write.table(as.data.frame(ego), 'Epi cell_identify\\EPISL(V2)DEG\\KEGG.txt', sep = '\t', row.names = FALSE, quote = FALSE)

# load package
library(ggplot2)
library(stringr)
# load data
go <- read.table(file = "go.select4.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
go <- go[c(5:1,9:6,14:10),]
go$Term <- factor(go$Description, levels = go$Description)
#go$Term <- go$Description

go$GeneRatio <- as.numeric(sapply(strsplit(go$GeneRatio,"/",fixed = T),"[",1))/as.numeric(sapply(strsplit(go$GeneRatio,"/",fixed = T),"[",2))
ggplot(go, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Term)) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(size = 10, colour = "black")) + 
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab(NULL) + facet_grid(ONTOLOGY~.,scales = "free") + 
  ggtitle("GO pathway enrichment") + scale_y_discrete(labels=function(x) str_wrap(x, width=40))
ggsave(filename = "go enrichment20211125V2全.pdf", width = 6,height = 4)



## Figure 2E ## 
library(ComplexHeatmap)
library(Seurat)
load("raw/seurat_epi.Rdata")
load("ZYJ/12/Epi_metadata20211125V2.Rdat")
seurat@meta.data <- Epi_metadata
mtx <- AverageExpression(seurat, group.by = "celltypes2")[[1]]
plotgenes <- c("GSTP1", "NQO1", "NQO2", "DUOX1", "DUOX2", "DUOXA1", "DUOXA2", "GPX2", "AQP5")
plotgenes%in%rownames(mtx)
all(plotgenes%in%rownames(mtx))
plot.data <- mtx[plotgenes,]

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}
plot.data <- standarize.fun(plot.data)

anncol<-data.frame(celltype = factor(colnames(plot.data),levels = colnames(plot.data)),
                   row.names = colnames(plot.data),
                   stringsAsFactors = F)
anncolors <- list()
anncolors[["celltype"]]<-c("Epi-SL" = "#074083",
                           # "Epi-3" = "#E10C0D",
                           # "Epi-1" = "#33A02C",
                           "Transit amplifying cells" = "#C9B1D5",
                           "Goblet cells" = "#FB9998",
                           "Intermediate" = "#FDBE6E",
                           "BEST4+ colonocytes" = "#A6CEE3",
                           "Tuft cells" = "#FF7F00",
                           "Enteroendocrine cells" = "#6A3D9A",
                           "Mature colonocytes" = "#B2DF8A") 
pdf("plot/Heatmap/Epi1202-4.pdf", width = 606/100, height = 494/100,)
pheatmap(plot.data,
         cluster_cols = F, cluster_rows = F,
         border_color = "white",
         # color = colorRampPalette(colors = c("#6353A3","white","#9A0038"))(100),
         color = colorRampPalette(colors = c("#75A3D0","white","red"))(100),
         annotation_col = anncol[colnames(plot.data),,drop=F],
         annotation_colors = anncolors,
         cellwidth = 12, cellheight = 12,
         treeheight_row = 20)
# gaps_row = 2)
invisible(dev.off())



## Figure 2G ## 
path <- "i:/genomicdata/external/ZYJ/11";setwd(path)

heatmap.BlWtRd <- c("#6699CC","white","#FF3C38")
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
heatmap.GrWtRd <- c("#2b2d42","#8d99ae","#edf2f4","#ef233c","#d90429")
heatmap.L.BlYlRd <- c("#4281a4","#9cafb7","#ead2ac","#e6b89c","#fe938c")
heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")
heatmap.fancy <- c("#10040A", "#2A0B35", "#4D155B", "#73215B", "#9C3558", "#C34D44", "#E07038", "#F2981C", "#F2CA51", "#FAF6A3")

expr <- openxlsx::read.xlsx("GSE76987_RightColonProcessed(1).xlsx")
plot.gene <- c("GSTP1", "NQO1", "NQO2", "DUOX1", "DUOX2", "DUOXA1", "DUOXA2", "AQP5")
plot.sample <- c(paste0("FPKM.CR-", seq(1:10)), 
                 paste0("FPKM.AP-", seq(1:10)), 
                 paste0("FPKM.SSA/P-", seq(1:21)))
plot.sample%in%colnames(expr)
plotdata <- expr[match(plot.gene, expr$gene), plot.sample]
rownames(plotdata) <- plot.gene
plotdata <- as.matrix(plotdata)
plotdata <- log(plotdata+1)
range(plotdata)
plotdata <- t(scale(t(plotdata)))
range(plotdata)
plotdata[plotdata<(-2)] <- -2
plotdata[plotdata>2] <- 2
range(plotdata)
plot(density(plotdata))
Group <- c(rep("NC", 10), rep("AP", 10), rep("SSL", 21))
anncol <- data.frame("Group" = factor(Group, levels = unique(Group)),
                     row.names = colnames(plotdata))
anncolors <- list()
anncolors[["Group"]] <- c("NC" = "#DFC27D",
                          "AP" = "#689E45",
                          "SSL" = "#5A8CA8")
# plot(x = as.factor(Group), y = plotdata[4,])

library(ComplexHeatmap)
pdf("GSE76987_FPKM_Heatmap.pdf", width = 1011/100, height = 443/100)
pheatmap(mat = plotdata, 
         cluster_rows = F, cluster_cols = F,
         border_color = "NA", gaps_col = c(10,20),
         # color = colorRampPalette(colors = c("#75A3D0","white","red"))(100),
         color = colorRampPalette(colors = heatmap.fancy)(100),
         cellwidth = 10, cellheight = 20,
         angle_col = "45",
         fontsize_col = 12, fontsize_row = 10,
         annotation_col = anncol, annotation_colors = anncolors,
         fontsize = 15,
         show_colnames = F,
         name = "scaled log(FPKM+1)")
invisible(dev.off())



## Figure 2H ## 
library(GSVA)
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

library(limma)
library(ggplot2)
library(ggpubr)
jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")

expr <- read.table("gse39582.expr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
pd <- read.table("easy_input_GSE39582.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
pd <- na.omit(pd)
pd <- data.frame(Samples = pd$GEOID,
                 Group = pd$MMR.Status,
                 stringsAsFactors = FALSE)
intersect(colnames(expr),pd$Samples)

go <- gmt2list("GOBP_RESPONSE_TO_OXIDATIVE_STRESS.gmt")
go.gsva <- gsva(expr = as.matrix(expr[,sinfo$GEOID]),
                go,
                method = "gsva")

sinfo$go <- as.numeric(go.gsva[1,])
kruskal.test(sinfo$go~sinfo$KRASBRAF)
boxplot(sinfo$go~sinfo$KRASBRAF)

sinfo$group <- ifelse(sinfo$KRASBRAF == 1,"BRAF-mut","WT")

ggplot(data = sinfo,aes(x = group, #分组列名
                        y = go, #连续变量列名
                        fill = group))+ #按分组填充颜色
  scale_fill_manual(values = jco[2:1]) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75),
              size = 0.8, color="white") +
  
  geom_point(shape = 21, size=2, # 点的性状和大小
             position = position_jitterdodge(), # 让点散开
             aes(color=group), alpha = 1) +
  scale_color_manual(values = jco[2:1]) + #用自定义颜色填充
  geom_boxplot(notch = TRUE, outlier.size = -1,
               color="black", lwd=0.8, alpha = 0.7) +
  theme_bw() + 
  ylab("RESPONSE_TO_OXIDATIVE_STRESS\nGSVA") +
  xlab("") +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  stat_compare_means(method = "wilcox.test")
ggsave("boxplot for go signature in gse39582 regarding BRAF.pdf",width = 3,height = 4)

