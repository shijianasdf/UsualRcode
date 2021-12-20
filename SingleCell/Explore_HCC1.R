library(ArchR)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(data.table)
library(dittoSeq)
library(Seurat)
library(stringr)
library(scCustomize)
library(GenomicRanges)
library(GenomicFeatures)
library(foreach)
library(doParallel)
library(data.table)
library(dplyr)
library(tibble)
library(future)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 50000 * 1024^2) #50G
library(patchwork)
library(cowplot)
library(ggsci)
library(Signac)
colors <- dittoColors()
zylcolor40<-c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF",
              "#7CFC00","#FFFF00","#808000","#FF00FF","#D2B48C","#7B68EE","#9400D3","#800080",
              "#A0522D","#FA8072","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD",
              "#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C",
              "#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
geneNum <- function(x){
  pos <- which(x > 0)
  return(length(pos))
}
cellNum <- function(x){
  pos <- which(x > 0)
  return(length(pos))
}
# 探索表达数据
{
  # the 10x hdf5 file contains both data types. cellranger细胞过滤后的数据
  hcc_1 <- Read10X_h5("/Volumes/HDD4/HCC_mulimodel/HCC_1/filtered_feature_bc_matrix.h5")
  df1 <- fread("/Volumes/HDD4/HCC_mulimodel/HCC_1/clusters.tsv")
  df2 <- fread("/Volumes/HDD4/HCC_mulimodel/HCC_1/donor_ids.tsv")
  df <- df1 %>% left_join(.,df2,by=c("barcode"="cell"))
  # extract RNA and ATAC data
  rna_counts <- hcc_1$`Gene Expression`
  atac_counts <- hcc_1$Peaks
  dim(rna_counts)
  dim(atac_counts)
  apply(rna_counts,2,geneNum) -> genesCell
  boxplot(genesCell)
  plot(density(genesCell))
  mean(genesCell)
  apply(rna_counts,1,cellNum) -> td
  td <- td[td != 0]
  boxplot(td)
  plot(density(td))
  mean(td)
  
  apply(atac_counts,2,geneNum) -> genesCell
  boxplot(genesCell)
  plot(density(genesCell))
  mean(genesCell)
  apply(atac_counts,1,cellNum) -> td
  td <- td[td != 0]
  boxplot(td)
  plot(density(td))
  mean(td)
  
  head(df1)
  common.anno <- getAnnotation(gtf.path = "/Users/biofly/project/scVelo/refdata-gex-GRCh38-2020-A/genes/genes.gtf")
  ##过滤双细胞,unssigned细胞
  rna_counts <- rna_counts[,df1$barcode] 
  pos <- which(df1$status == "singlet")
  rna_counts <- rna_counts[,pos]
  df1 <- df1[pos,]
  
  sce <- FastSeurat_1(count=rna_counts,
                    dir.name=NULL,
                    obj=NULL,#min.cells=3,
                    species=c("human","mouse")[1],
                    min.features=0,
                    max.features=300000,
                    percent.mt.num=100,
                    plot=F,
                    pcSelect=50,
                    project="HCC_1",
                    nfeatures=2000,
                    vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                    all.scale = F,
                    npcs = 50,
                    resolution = 0.5,
                    harmony=F,
                    doublet=F,
                    perplexity = 30,
                    isMarkers=T,
                    cellCycle=F,
                    rmOtherGene=T,
                    features=NULL,
                    filepath=NULL,
                    outdir="Results",
                    names="love")
  sce$v
  bardata <- sce$sce@meta.data %>% 
    dplyr::group_by(orig.ident) %>% 
    dplyr::summarise(num=n())
  p0 <- ggplot(data=bardata,mapping=aes(x=orig.ident,y=num,fill=orig.ident)) +
    geom_bar(stat = "identity") +
    #geom_text_repel() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
          axis.text=element_text(face = "bold",size = 14),
          legend.text=element_text(face = "bold",size = 14),
          legend.title = element_text(face = "bold",size = 14),
          axis.title = element_text(face = "bold",size = 14))
  VlnPlot(sce$sce,features = c('nCount_RNA','nFeature_RNA','percent.mt',"percent.ribo"),group.by = "orig.ident",pt.size=0)
  head(sce$sce@meta.data)
  sce$sce@meta.data$orig.ident <- df1$assignment
  match(rownames(sce$sce@meta.data), df1$barcode)
  p2 <- dittoDimPlot(sce$sce, reduction.use  = "umap",var="seurat_clusters",size = 1,do.label = TRUE, do.ellipse = TRUE,legend.size = 9,shape.legend.size=9,labels.size=7,do.raster=T,raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 16))
  p3 <- dittoDimPlot(sce$sce, reduction.use = "umap",var="orig.ident",size = 1,do.label = TRUE, do.ellipse = TRUE,legend.size = 9,shape.legend.size=9,labels.size=7,do.raster=T,raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 16))
  p2 | p3 
  p4 <- VlnPlot(sce$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  cl_markers <- sce$sce.markers
  cl_markers <- cl_markers$gene[convert(cl_markers$gene,fromtype = "SYMBOL",totype="gene_type",db=common.anno)=="protein_coding"]
  top10_cl_markers <- cl_markers %>% dplyr::group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
  dh <- DoHeatmap(sce$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  #EPCAM <- c("EPCAM","KRT19","KRT7") #上皮细胞
  #cholangiocytes <- c("FYXD2","TM4SF4", "ANXA4") #胆管细胞
  progenitor.cell <- c("KRT19")
  hepatocytes  <- c("ALB","KRT8") # 肝脏细胞
  fibroblasts <- c("ACTA2" , "THY1") #成纤维细胞
  endothelial <- c("PECAM1","ENG","VWF")
  All.Immune <- "PTPRC" #免疫细胞
  B.cell <- c("MZB1","CD79A","CD79B","MS4A1")
  T.cell <- c("CD3D","CD3E","CD8A","CD4","IL7R")
  myeloid <- c("LYZ","CD14")
  Dendritic.cell <- c("CLEC9A","CD1C","THBD","ITGAX","CD1C")
  Macrophage <- c("CD68","CD14","HLA-DRA","HLA-DRB1","APOE","MMP9","CTSK")
  Mast <- c("KIT","TPSAB1","TPSB2","ENPP3")
  Monocyte <- c("LYZ","CD68","CD163","ITGAX","S100A8","S100A9")
  NK <- c("GNLY","NKG7")
  pDC <- c("LILRA4","UGCG","CLEC4C")
  plasma <- c("SDC1","CD38","SLAMF7","DERL3","MZB1","IGHA1","IGHA2")
  proliferative <- c("MKI67","PCNA","TK1","TOP2A","TUBB","TYMS")
  FeaturePlot_scCustom(sce$sce, features = progenitor.cell, order = F)
  FeaturePlot_scCustom(sce$sce, features = EPCAM, order = F)
  FeaturePlot_scCustom(sce$sce, features = cholangiocytes, order = F)
  FeaturePlot_scCustom(sce$sce, features = hepatocytes,order = F)
  FeaturePlot_scCustom(sce$sce, features = fibroblasts,order = F)
  FeaturePlot_scCustom(sce$sce, features = endothelial,order = F)
  FeaturePlot_scCustom(sce$sce, features = All.Immune,order = F)
  FeaturePlot_scCustom(sce$sce, features = B.cell,order = F)
  FeaturePlot_scCustom(sce$sce, features = plasma,order = F)
  FeaturePlot_scCustom(sce$sce, features = T.cell,order = F)
  FeaturePlot_scCustom(sce$sce, features = myeloid,order = F)
  FeaturePlot_scCustom(sce$sce, features = NK,order = F)
  FeaturePlot_scCustom(sce$sce, features = Dendritic.cell,order = F)
  FeaturePlot_scCustom(sce$sce, features = Macrophage,order = F)
  FeaturePlot_scCustom(sce$sce, features = Mast,order = F)
  FeaturePlot_scCustom(sce$sce, features = Monocyte,order = F)
  FeaturePlot_scCustom(sce$sce, features = pDC,order = F)
  #infercnv pipeline
  normal <- c(7,10,11,13)
  non_normal <- c(0,1,2,3,4,5,6,8,9,12,14,15)
  sce$sce@meta.data$celltypes <- NA
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==7] ="immue cell_7"
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==10] = "immue cell_10"
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==11] = "fibroblasts_11"
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==13] = "endothelial_13"
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==0] = "hepatocytes_0"
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==1] = "hepatocytes_1"
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==2] = "hepatocytes_2"
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==3] = "hepatocytes_3"
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==4] = "hepatocytes_4"
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==5] = "hepatocytes_5"
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==6] = "hepatocytes_6"
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==8] = "hepatocytes_8"
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==9] = "hepatocytes_9"
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==12] = "hepatocytes_12"
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==14] = "hepatocytes_14"
  sce$sce@meta.data$celltypes[sce$sce@meta.data$seurat_clusters==15] = "hepatocytes_15"
  res <- FastInferCNV(obj=sce$sce,
                      counts=NULL,
                      annotations_file=NULL,
                      cellType="celltypes",
                      gene_order_file=NULL,
                      ref_group_names=c("immue cell_7","immue cell_10","fibroblasts_11","endothelial_13"),  
                      cutoff=c(0.1,1)[1], 
                      analysis_mode=c("samples", "subclusters", "cells")[1],
                      out_dir="HCC_1_inferCNV",
                      cluster=T,
                      HMM=F,
                      denoise=T,
                      anno.db="/Users/biofly/project/scVelo/refdata-gex-GRCh38-2020-A/genes/genes.gtf")
  
  sce <- FastSeurat(count=rna_counts,
                    dir.name=NULL,
                    obj=NULL,#min.cells=3,
                    species=c("human","mouse")[1],
                    min.features=200,
                    max.features=10000,
                    percent.mt.num=50,
                    plot=F,
                    pcSelect=50,
                    project="HCC_1",
                    nfeatures=2000,
                    vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                    all.scale = F,
                    npcs = 50,
                    resolution = 0.5,
                    harmony=F,
                    doublet=T,
                    perplexity = 30,
                    isMarkers=T,
                    cellCycle=T,
                    rmOtherGene=F,
                    features=NULL,
                    filepath=NULL,
                    outdir="Results",
                    names="love")
  sce$v
  p2 <- dittoDimPlot(sce$sce, reduction.use  = "umap",var="seurat_clusters",size = 1,do.label = TRUE, do.ellipse = TRUE,legend.size = 9,shape.legend.size=9,labels.size=7,do.raster=T,raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 16))
  p2 
  p4 <- VlnPlot(sce$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  EPCAM <- c("EPCAM","KRT19","KRT7") #上皮细胞
  cholangiocytes <- c("FYXD2","TM4SF4", "ANXA4") #胆管细胞
  hepatocytes  <- c("APOC3","FABP1","APOA1","ALB","KRT8") # 肝脏细胞
  fibroblasts <- c("ACTA2" , "THY1", "COL1A2") #成纤维细胞
  endothelial <- c("PECAM1","ENG","VWF")
  All.Immune <- "PTPRC" #免疫细胞
  B.cell <- c("MZB1","CD79A","CD79B","MS4A1")
  T.cell <- c("CD3D","CD3E","CD8A","CD4","IL7R")
  myeloid <- c("LYZ","CD14")
  Dendritic.cell <- c("CLEC9A","CD1C","THBD","ITGAX","CD1C")
  Macrophage <- c("CD68","CD14","HLA-DRA","HLA-DRB1","APOE","MMP9","CTSK")
  Mast <- c("KIT","TPSAB1","TPSB2","ENPP3")
  Monocyte <- c("LYZ","CD68","CD163","ITGAX","S100A8","S100A9")
  NK <- c("GNLY","NKG7")
  pDC <- c("LILRA4","UGCG","CLEC4C")
  plasma <- c("SDC1","CD38","SLAMF7","DERL3","MZB1","IGHA1","IGHA2")
  proliferative <- c("MKI67","PCNA","TK1","TOP2A","TUBB","TYMS")
  FeaturePlot_scCustom(sce$sce, features = EPCAM, order = F)
  FeaturePlot_scCustom(sce$sce, features = cholangiocytes, order = F)
  FeaturePlot_scCustom(sce$sce, features = hepatocytes,order = F)
  FeaturePlot_scCustom(sce$sce, features = fibroblasts,order = F)
  FeaturePlot_scCustom(sce$sce, features = endothelial,order = F)
  FeaturePlot_scCustom(sce$sce, features = All.Immune,order = F)
  FeaturePlot_scCustom(sce$sce, features = B.cell,order = F)
  FeaturePlot_scCustom(sce$sce, features = plasma,order = F)
  FeaturePlot_scCustom(sce$sce, features = T.cell,order = F)
  FeaturePlot_scCustom(sce$sce, features = proliferative,order = F)
  FeaturePlot_scCustom(sce$sce, features = NK,order = F)
  FeaturePlot_scCustom(sce$sce, features = Dendritic.cell,order = F)
  FeaturePlot_scCustom(sce$sce, features = Macrophage,order = F)
  FeaturePlot_scCustom(sce$sce, features = Mast,order = F)
  FeaturePlot_scCustom(sce$sce, features = Monocyte,order = F)
  FeaturePlot_scCustom(sce$sce, features = pDC,order = F)
  cl_markers <- sce$sce.markers
  top10_cl_markers <- cl_markers %>% dplyr::group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
  dh <- DoHeatmap(sce$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  
  rna_counts <- sce$sce@assays$RNA@counts
  dim(rna_counts)
  apply(rna_counts,2,geneNum) -> genesCell
  boxplot(genesCell)
  plot(density(genesCell))
  mean(genesCell)
  apply(rna_counts,1,cellNum) -> td
  td <- td[td != 0]
  boxplot(td)
  plot(density(td))
  mean(td)
  
  apply(atac_counts,2,geneNum) -> genesCell
  boxplot(genesCell)
  plot(density(genesCell))
  mean(genesCell)
  apply(atac_counts,1,cellNum) -> td
  td <- td[td != 0]
  boxplot(td)
  plot(density(td))
  mean(td)
}

#探索atac数据
{
  addArchRThreads(threads = 16)
  addArchRGenome("hg38")
  genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Hsapiens.UCSC.hg38)
  geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, OrgDb = org.Hs.eg.db)
  addArchRChrPrefix(chrPrefix = FALSE)
  frag.file <- "/Volumes/HDD4/HCC_mulimodel/HCC_1/atac_fragments.tsv.gz"
  #Create arrow files
  ArrowFiles <- createArrowFiles(inputFiles = frag.file,sampleNames = "HCC_1",minTSS = 0, minFrags = 0, maxFrags=1e+10000, excludeChr = NULL,addTileMat = TRUE, addGeneScoreMat = TRUE,force = T)
  #Add doublet scores
  doubScores <- addDoubletScores(input = ArrowFiles, k = 10, knnMethod = "UMAP", LSIMethod = 1)
  #Create arrow project
  proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "HCC_atac", copyArrows = TRUE)
  
  #过滤ArchR对象的细胞(保留表达矩阵中的细胞)
  #fragments <- fread("/Volumes/HDD4/HCC_mulimodel/HCC_1/atac_fragments.tsv.gz")
  #length(unique(fragments$V4))
  df1 <- fread("/Volumes/HDD4/HCC_mulimodel/HCC_1/clusters.tsv")
  pos <- which(df1$status == "singlet")
  df1 <- df1[pos,]
  cellid <- str_split(proj$cellNames,"#",simplify = T)[,2]
  head(df1);dim(df1)
  pos <- match(df1$barcode,cellid)
  proj <- proj[proj$cellNames[pos],]
  paste0("Memory Size = ", round(object.size(proj) / 10^6, 3), "MB")
  getAvailableMatrices(proj)
  # remove doublets
  #proj <- filterDoublets(proj, filterRatio = 1.5)
  saveArchRProject(ArchRProj = proj)
  #线性降维
  proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2, clusterParams = list(resolution = c(0.2),sampleCells = 10000, n.start = 10), varFeatures = 25000, dimsToUse = 1:30)
  #Clustering using Seurat’s FindClusters() function
  proj <- addClusters(input = proj,reducedDims = "IterativeLSI",method = "Seurat",name = "Clusters",resolution = 0.8)
  #非线性降维可视化
  proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine")
  proj <- addTSNE(ArchRProj = proj, reducedDims = "IterativeLSI", name = "TSNE", perplexity = 30)
  
  head(proj@cellColData)
  cellid <- str_split(proj$cellNames,"#",simplify = T)[,2]
  match(cellid,df1$barcode)
  proj@cellColData$Sample <- df1$assignment
  proj@cellColData$Sample <- rep("HCC_1",length(proj@cellColData$Sample))
  proj@cellColData$sampleID <- df1$assignment
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "sampleID", embedding = "UMAP",size = 0.1,baseSize = 15)
  p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP",size = 0.1,baseSize = 15) 
  ggAlignPlots(p1, p2, type = "h")
  plotFragmentSizes(ArchRProj = proj)
  plotTSSEnrichment(ArchRProj = proj)
  df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
  ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
  ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
  
  #保存ArchR对象
  subsetArchRProject(ArchRProj=proj, cells = getCellNames(proj), outputDirectory = "ArchRSubset")
  #重新加载对象
  proj <- loadArchRProject(path="ArchRSubset")
  #识别marker基因
  markersGS <- getMarkerFeatures(ArchRProj = proj, useMatrix = "GeneScoreMatrix", groupBy = "Clusters",bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
  markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
  markerList$C6
  heatmapGS <- plotMarkerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1",transpose = TRUE)
  ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  proj <- addImputeWeights(proj)
  #细胞类型注释
  progenitor.cell <- c("KRT19")
  hepatocytes  <- c("ALB","KRT8") # 肝脏细胞
  fibroblasts <- c("ACTA2" , "THY1") #成纤维细胞
  endothelial <- c("PECAM1","VWF")
  All.Immune <- "PTPRC" #免疫细胞
  B.cell <- c("MZB1","CD79A","CD79B","MS4A1")
  T.cell <- c("CD3D","CD3E","CD8A","CD4","IL7R")
  NK <- c("GNLY","NKG7")
  myeloid <- c("LYZ","CD14")
  Dendritic.cell <- c("CLEC9A","CD1C","THBD","ITGAX","CD1C")
  Macrophage <- c("CD68","CD14","HLA-DRA","HLA-DRB1","APOE","MMP9","CTSK")
  Mast <- c("KIT","TPSAB1","TPSB2","ENPP3")
  Monocyte <- c("LYZ","CD68","CD163","ITGAX","S100A8","S100A9")
  pDC <- c("LILRA4","UGCG","CLEC4C")
  plasma <- c("SDC1","CD38","SLAMF7","DERL3","MZB1","IGHA1","IGHA2")
  proliferative <- c("MKI67","PCNA","TK1","TOP2A","TUBB","TYMS")
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = progenitor.cell, embedding = "UMAP",imputeWeights = getImputeWeights(proj),size = 0.1,baseSize = 15)
  p3 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = hepatocytes, embedding = "UMAP",imputeWeights = getImputeWeights(proj),size = 0.1,baseSize = 15)
  p4 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = fibroblasts, embedding = "UMAP",imputeWeights = getImputeWeights(proj),size = 0.1,baseSize = 15)
  p5 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = endothelial, embedding = "UMAP",imputeWeights = getImputeWeights(proj),size = 0.1,baseSize = 15)
  p6 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = All.Immune, embedding = "UMAP",imputeWeights = getImputeWeights(proj),size = 0.1,baseSize = 15)
  p7 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = B.cell, embedding = "UMAP",imputeWeights = getImputeWeights(proj),size = 0.1,baseSize = 15)
  p8 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = T.cell, embedding = "UMAP",imputeWeights = getImputeWeights(proj),size = 0.1,baseSize = 15)
  p9 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = NK, embedding = "UMAP",imputeWeights = getImputeWeights(proj),size = 0.1,baseSize = 15)
  p10 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = Macrophage, embedding = "UMAP",imputeWeights = getImputeWeights(proj),size = 0.1,baseSize = 15)
  p11 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = Monocyte, embedding = "UMAP",imputeWeights = getImputeWeights(proj),size = 0.1,baseSize = 15)
  p12 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = Dendritic.cell, embedding = "UMAP",imputeWeights = getImputeWeights(proj),size = 0.1,baseSize = 15)
  p13 <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = myeloid, embedding = "UMAP",imputeWeights = getImputeWeights(proj),size = 0.1,baseSize = 15)
  pp <- lapply(p12, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 15) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
      )
  })
  do.call(cowplot::plot_grid, c(list(ncol = 3),pp))
  
}

#WNN 数据聚类分群
{
  hcc_1 <- Read10X_h5("/Volumes/HDD4/HCC_mulimodel/HCC_1/filtered_feature_bc_matrix.h5")
  # extract RNA and ATAC data
  rna_counts <- hcc_1$`Gene Expression`
  atac_counts <- hcc_1$Peaks
  # Create Seurat object
  pbmc <- CreateSeuratObject(counts = rna_counts)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  # Now add in the ATAC-seq data
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "hg38"
  frag.file <- "/Volumes/HDD4/HCC_mulimodel/HCC_1/atac_fragments.tsv.gz"
  chrom_assay <- CreateChromatinAssay(counts = atac_counts,sep = c(":", "-"),genome = 'hg38',fragments = frag.file,min.cells = 10,annotation = annotations)
  pbmc[["ATAC"]] <- chrom_assay
  #细胞质量控制
  VlnPlot(pbmc, features = c("nCount_ATAC", "nFeature_RNA", "nCount_RNA","nFeature_ATAC","percent.mt"), ncol = 3,
          log = TRUE, pt.size = 0) + NoLegend()
  pbmc <- subset(x = pbmc,subset = nFeature_RNA < 10000 & nFeature_RNA > 200 & percent.mt < 50)
  # RNA analysis
  DefaultAssay(pbmc) <- "RNA"
  pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% FindNeighbors() %>% FindClusters() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  # ATAC analysiss
  # We exclude the first dimension as this is typically correlated with sequencing depth
  DefaultAssay(pbmc) <- "ATAC"
  pbmc <- RunTFIDF(pbmc)
  pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
  pbmc <- RunSVD(pbmc)
  pbmc <- FindNeighbors(pbmc, reduction = 'lsi', dims = 2:50)
  pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.8)
  pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  #我们计算 WNN 图，表示 RNA 和 ATAC-seq 模式的加权组合。我们将此图用于 UMAP 可视化和聚类
  pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
  pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
  
  #p3 <- dittoDimPlot(pbmc, reduction.use  = "wnn.umap",var="seurat_clusters",size = 1,do.label = TRUE, do.ellipse = TRUE,legend.size = 9,shape.legend.size=9,labels.size=7,do.raster=T,raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 16))
  p4 <- dittoDimPlot(pbmc, reduction.use  = "wnn.umap",var="wsnn_res.0.8",size = 1,do.label = TRUE, do.ellipse = TRUE,legend.size = 9,shape.legend.size=9,labels.size=7,do.raster=T,raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 16))
  #p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  p4
  
  EPCAM <- c("EPCAM","KRT19","KRT7") #上皮细胞
  cholangiocytes <- c("FYXD2","TM4SF4", "ANXA4") #胆管细胞
  hepatocytes  <- c("APOC3","FABP1","APOA1") # 肝脏细胞
  fibroblasts <- c("ACTA2" ,"COL1A2") #成纤维细胞
  endothelial <- c("ENG","VWF")
  All.Immune <- "PTPRC" #免疫细胞
  B.cell <- c("BANK1","CD19","CD79A","CD79B","MS4A1")
  T.cell <- c("CD3D","CD3E","CD8A","CD4","CD2")
  Dendritic.cell <- c("CLEC9A","CD1C","THBD","ITGAX","CD1C")
  Macrophage <- c("CD68","CD14","HLA-DRA","HLA-DRB1","APOE","MMP9","CTSK")
  Mast <- c("KIT","TPSAB1","TPSB2","ENPP3")
  Monocyte <- c("LYZ","CD68","CD163","ITGAX","S100A8","S100A9")
  NK <- c("FCGR3A","FCGR3B","NCAM1","PVRIG","TIGIT","CD7", "FGFBP2", "KLRF1")
  pDC <- c("LILRA4","UGCG","CLEC4C")
  plasma <- c("SDC1","CD38","SLAMF7","DERL3","MZB1","IGHA1","IGHA2")
  proliferative <- c("MKI67","PCNA","TK1","TOP2A","TUBB","TYMS")
  DefaultAssay(pbmc) <- "RNA"
  FeaturePlot_scCustom(pbmc, features = EPCAM, reduction = 'wnn.umap',order = F)
  FeaturePlot_scCustom(pbmc, features = cholangiocytes, reduction = 'wnn.umap',order = F)
  FeaturePlot_scCustom(pbmc, features = hepatocytes,reduction = 'wnn.umap',order = F)
  FeaturePlot_scCustom(pbmc, features = fibroblasts,reduction = 'wnn.umap',order = F)
  FeaturePlot_scCustom(pbmc, features = endothelial,reduction = 'wnn.umap',order = F)
  FeaturePlot_scCustom(pbmc, features = All.Immune,reduction = 'wnn.umap',order = F)
  FeaturePlot_scCustom(pbmc, features = B.cell,reduction = 'wnn.umap',order = F)
  FeaturePlot_scCustom(pbmc, features = plasma,reduction = 'wnn.umap',order = F)
  FeaturePlot_scCustom(pbmc, features = T.cell,reduction = 'wnn.umap',order = F)
  FeaturePlot_scCustom(pbmc, features = proliferative,reduction = 'wnn.umap',order = F)
  FeaturePlot_scCustom(pbmc, features = NK,reduction = 'wnn.umap',order = F)
  FeaturePlot_scCustom(pbmc, features = Dendritic.cell,reduction = 'wnn.umap',order = F)
  FeaturePlot_scCustom(pbmc, features = Macrophage,reduction = 'wnn.umap',order = F)
  FeaturePlot_scCustom(pbmc, features = Mast,reduction = 'wnn.umap',order = F)
  FeaturePlot_scCustom(pbmc, features = Monocyte,reduction = 'wnn.umap',order = F)
  FeaturePlot_scCustom(pbmc, features = pDC,reduction = 'wnn.umap',order = F)
  
  #Seurat计算marker
  DefaultAssay(pbmc) <- "RNA"
  markers1 <- FindAllMarkers(pbmc)
  cl_markers <- markers1
  top10_cl_markers1 <- cl_markers %>% dplyr::group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group=T) %>% top_n(n = 10, wt = avg_log2FC) 
  dh <- DoHeatmap(pbmc, features = top10_cl_markers1$gene,slot = "data") + NoLegend()
  print(dh)
  #presto计算marker
  library(presto)
  markers <- wilcoxauc(pbmc)
  markers %>%
    group_by(group) %>%
    arrange(desc(logFC), .by_group=T) %>%
    top_n(n = 10, wt = logFC) -> top10_cl_markers
  dh <- DoHeatmap(pbmc, features = top10_cl_markers$feature, slot = "data") + NoLegend()
  print(dh)
}
