library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "data/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")
str(pbmc.data)
dense.size <- object.size(x = as.matrix(x = pbmc.data))
dense.size
sparse.size <- object.size(x = pbmc.data)
sparse.size
dense.size/sparse.size
# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           min.cells = 3, 
                           min.features  = 200, 
                           project = "10X_PBMC")
str(pbmc)
# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
mito.genes <- grep(pattern = "^MT-", x = rownames(pbmc@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@assays[["RNA"]][mito.genes, ])/Matrix::colSums(pbmc@assays[["RNA"]])
# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
#Seurat v2 function, but shows compatibility in Seurat v3
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
#pbmc$percent.mito <- percent.mito
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
par(mfrow = c(1, 2))
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# We filter out cells that have unique gene counts (nFeature_RNA) over 2,500 or less than
# 200 Note that > and < are used to define a'gate'.  
#-Inf and Inf should be used if you don't want a lower or upper threshold.
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito >  -Inf & percent.mito < 0.05 )

#After removing unwanted cells from the dataset, the next step is to normalize the data.
{
  pbmc@assays$RNA@data
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000) # 对count矩阵进行lognomalise
  pbmc@assays$RNA@data
}

#Seurat calculates highly variable genes and focuses on these for downstream analysis.
{
  pbmc <- FindVariableFeatures(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
  #To view the output of the FindVariableFeatures output we use this function. 
  head(x = HVFInfo(object = pbmc))
  head(VariableFeatures(pbmc), 10)
  pbmc@assays$RNA@var.features
  pbmc@assays$RNA@meta.features
}

#Scaling the data and removing unwanted sources of variation
{
  #Your single cell dataset likely contains ‘uninteresting’ sources of variation.
  #This could include not only technical noise, but batch effects, or even biological sources of variation (cell cycle stage).
  #To mitigate the effect of these signals, Seurat constructs linear models to predict gene expression based on user-defined variables. 
  #The scaled z-scored residuals of these models are stored in the scale.data slot, 
  #and are used for dimensionality reduction and clustering.
  pbmc@assays$RNA@scale.data[1:4,1:4]
  pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nCounts_RNA", "percent.mito"))
  pbmc@assays$RNA@scale.data[1:4,1:4]
}

#Perform linear dimensional reduction
{
  pbmc <- RunPCA(object = pbmc,  npcs = 30, verbose = FALSE)
  pbmc@reductions$pca@feature.loadings
  pbmc@reductions$pca@cell.embeddings
  # Examine and visualize PCA results a few different ways
  DimPlot(object = pbmc, reduction = "pca")
  # Dimensional reduction plot, with cells colored by a quantitative feature
  FeaturePlot(object = pbmc, features = "MS4A1")
  # Scatter plot across single cells, replaces GenePlot
  FeatureScatter(object = pbmc, feature1 = "MS4A1", feature2 = "PC_1")
  FeatureScatter(object = pbmc, feature1 = "MS4A1", feature2 = "CD3D")
  # Violin and Ridge plots
  VlnPlot(object = pbmc, features = c("LYZ", "CCL5", "IL32"))
  RidgePlot(object = pbmc, feature = c("LYZ", "CCL5", "IL32"))
  # Heatmaps 
  #Both cells and genes are ordered according to their PCA scores. 奇怪这个图代表什么
  DimHeatmap(object = pbmc, reduction = "pca", cells = 200, balanced = TRUE)
  
}

#Determine statistically significant principal components
{
  #In Macosko et al, we implemented a resampling test inspired by the jackStraw procedure. 
  #We randomly permute a subset of the data (1% by default) and rerun PCA, 
  #constructing a ‘null distribution’ of gene scores, and repeat this procedure. 
  #We identify ‘significant’ PCs as those who have a strong enrichment of low p-value genes.
  pbmc <- JackStraw(object = pbmc,reduction = "pca", dims = 20, num.replicate = 100,  prop.freq = 0.1, verbose = FALSE)
  pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20, reduction = "pca")
  JackStrawPlot(object = pbmc, dims = 1:20, reduction = "pca")
  ElbowPlot(object = pbmc)
  pbmc@reductions$pca@jackstraw$overall.p.values
}

#Cluster the cells
{
  pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:20)
  pbmc@neighbors
  pbmc <- FindClusters(pbmc, resolution = 0.5, algorithm = 1)
  pbmc@graphs$RNA_snn
  pbmc@graphs$RNA_nn
  pbmc@active.ident
  Idents(object = pbmc)
}

#Run Non-linear dimensional reduction (tSNE)
{
  pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE,perplexity=30)
  pbmc@reductions$tsne@cell.embeddings
  Idents(object = pbmc)
  DimPlot(object = pbmc, reduction = "tsne")
}

#Run UMAP
{
  pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:20)
  pbmc@reductions$umap@cell.embeddings
  DimPlot(pbmc, reduction = "umap")
  pbmc$tsne[[1:nrow(pbmc$tsne)]]
  pbmc$pca[[1:nrow(pbmc$pca)]]
  pbmc$umap[[1:nrow(pbmc$umap)]]
  Idents(object = pbmc)
  pbmc$seurat_clusters
  saveRDS(pbmc, file = "data/pbmc_tutorial.rds")
}

#Finding differentially expressed genes (cluster biomarkers)
{
  #Seurat can help you find markers that define clusters via differential expression
  # find all markers of cluster 1
  cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
  print(x = head(x = cluster1.markers, n = 5))
  # find all markers distinguishing cluster 2 from clusters 0 and 3
  cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 2, ident.2 = c(0, 3), min.pct = 0.25)
  print(x = head(x = cluster5.markers, n = 5))
  # find markers for every cluster compared to all remaining cells, report
  # only the positive ones
  pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
  cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)
  VlnPlot(object = pbmc, features =c("NKG7", "PF4"))
  FeaturePlot(object = pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols = c("grey", "blue"), reduction = "tsne")
  top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
  # setting slim.col.label to TRUE will print just the cluster IDS instead of
  # every cell name
  DoHeatmap(object = pbmc, features = top10$gene, label = TRUE)
}

#Assigning cell type identity to clusters
{
  new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                       "NK", "DC", "Platelet")
  names(new.cluster.ids) <- levels(pbmc)
  pbmc <- RenameIdents(pbmc, new.cluster.ids)
  DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
}

#Further subdivisions within cell types
{
  # First lets stash our identities for later
  pbmc <- StashIdent(object = pbmc, save.name = "ClusterNames_0.6")
  # Note that if you set save.snn=T above, you don't need to recalculate the
  # SNN, and can simply put: pbmc <- FindClusters(pbmc,resolution = 0.8)
  pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, resolution = 0.8, print.output = FALSE)
  # Demonstration of how to plot two tSNE plots side by side, and how to color
  # points based on different criteria
  plot1 <- DimPlot(object = pbmc, reduction = "tsne")
  plot2 <- DimPlot(object = pbmc, reduction = "tsne",group.by = "ClusterNames_0.6",)
  plot_grid(plot1, plot2)
  # Find discriminating markers
  tcell.markers <- FindMarkers(object = pbmc, ident.1 = 0, ident.2 = 1)
  # Most of the markers tend to be expressed in C1 (i.e. S100A4). However, we
  # can see that CCR7 is upregulated in C0, strongly indicating that we can
  # differentiate memory from naive CD4 cells.  cols.use demarcates the color
  # palette from low to high expression
  FeaturePlot(object = pbmc, features = c("S100A4", "CCR7"), cols = c("green", "blue"))
  pbmc <- SetIdent(object = pbmc, value = "ClusterNames_0.6")
  saveRDS(pbmc, file = "data/pbmc3k_final.rds")
}
