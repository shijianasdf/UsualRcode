#' ICC数据质控
#' @author shijian
#' 数据地址 ./data/peng_ICC
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)
library(dittoSeq)
library(scCustomize)
library(ggrepel)
library(reshape2)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
zylcolor40<-c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF",
              "#7CFC00","#FFFF00","#808000","#FF00FF","#D2B48C","#7B68EE","#9400D3","#800080",
              "#A0522D","#FA8072","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD",
              "#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C",
              "#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
scales::show_col(zylcolor40)
#天津肿瘤P1_653303样本数据质控
{
  sce <- FastSeurat(count=NULL,
                    dir.name="/Volumes/HDD3/tianjin_ICC/P1_653303/mRNA/outs/filtered_feature_bc_matrix",
                    obj=NULL,
                    species=c("human","mouse")[1],
                    min.features=200,
                    max.features=7500,
                    percent.mt.num=50,
                    plot=F,
                    pcSelect=30,
                    project="P1_653303",
                    nfeatures=2000,
                    vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                    all.scale = F,
                    npcs = 50,
                    resolution = 0.5,
                    harmony=F,
                    doublet=T,
                    isMarkers=T,
                    rmOtherGene=T,
                    perplexity = 30,
                    cellCycle=T,
                    features=NULL,
                    filepath=NULL,
                    outdir="Results",
                    names="love")
  save(sce,file="./Results/P1_653303_seurat.rda")
  sce$initial_sce #8188个细胞
  sce$double_sce #7900个细胞
  sce$sce #7308个细胞 592个双细胞
  print(sce$v) #min.features=200,max.features=7500,percent.mt.num=50
  print(sce$e) #pcSelect=30
  
  sce$sce #3203个细胞 2963个单细胞 240个双细胞
  print(sce$v) #min.features=200,max.features=4000,percent.mt.num=50
  print(sce$e) #pcSelect=30
  
  p1 <- sce$v
  p2 <- sce$e
  p3 <- DimPlot(sce$sce, reduction = "umap",label.size = 4,repel=T,label=T)
  p4 <- VlnPlot(sce$sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), idents=NULL,ncol = 3,pt.size=0)
  p2+p1+p3+p4+plot_layout(nrow=2)
  
  library(dplyr)
  #sce$sce.markers <- FindAllMarkers(object = sce$sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
  cl_markers <- sce$sce.markers
  top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  dh <- DoHeatmap(sce$sce, features = top10_cl_markers$gene) + NoLegend()
  
}

# 5个样本合并
{
  filenames <- list.files("./data/peng_ICC")
  files <- file.path(list.files("./data/peng_ICC",full.names = T),"filtered_feature_bc_matrix")
  sceList <- list()
  for(i in 1:length(files)){
    sceTemp <- FastSeurat(count=NULL,
                      dir.name=files[i],
                      obj=NULL,#min.cells=3,
                      min.features=200,
                      max.features=7500,
                      percent.mt.num=50,
                      plot=F,
                      pcSelect=30,
                      project=filenames[i],
                      nfeatures=2000,
                      vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                      all.scale = F,
                      npcs = 50,
                      resolution = 0.5,
                      harmony=F,
                      doublet=T,
                      isMarkers=F,
                      perplexity = 30,
                      cellCycle=F,
                      features=NULL,
                      filepath=NULL,
                      outdir="Results",
                      names="love")
    sceList[[i]] <- sceTemp
  }
  names(sceList) <- filenames
  sce.list <- lapply(sceList,function(x){
    x$sce
  })
  sce.list <- c(sce.list,P1_653303=sce$sce)
  sce_merge <- FastMergeSeurat(objList=sce.list,
                  object.names=c(filenames,"P1_653303"),
                  project.name.default="ICC")
  sce_merge_res <- FastSeurat(count=NULL,
                               dir.name=NULL,
                               obj=sce_merge,#min.cells=3,
                               min.features=200,
                               max.features=7500,
                               percent.mt.num=50,
                               plot=F,
                               pcSelect=30,
                               project="ICC",
                               nfeatures=2000,
                               vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                               all.scale = F,
                               npcs = 50,
                               resolution = 0.5,
                               harmony=F,
                               doublet=F,
                               perplexity = 30,
                               isMarkers=F,
                               cellCycle=T,
                               features=NULL,
                               filepath=NULL,
                               outdir="Results",
                               names="love")
  save(sce_merge_res,file="./Results/P1_653303_sce_merge_noHamony_res.rda")
  sce_merge_res$sce@reductions$harmony@cell.embeddings
  sce_merge_res$sce@meta.data
  p1 <- sce_merge_res$v
  
  #p4 <- VlnPlot(sce_merge_res$sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), idents=NULL,ncol = 3,pt.size=0)
  head(sce_merge_res$sce@meta.data)
  bardata <- sce_merge_res$sce@meta.data %>% 
                           group_by(orig.ident,Doublet) %>% 
                           summarise(num=n())
  p1 <- ggplot(data=bardata,mapping=aes(x=orig.ident,y=num,fill=Doublet)) +
            geom_bar(stat = "identity") +
            #geom_text_repel() +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  p2 <- sce_merge_res$e
  p3 <- DimPlot(sce_merge_res$sce, reduction = "umap",group.by="seurat_clusters",label.size = 4,repel=T,label=T)
  p4 <- DimPlot(sce_merge_res$sce, reduction = "umap",group.by="orig.ident",label.size = 4,repel=T,label=T)
  p5 <- VlnPlot(sce_merge_res$sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), idents=NULL,ncol = 3,pt.size=0)
  (p1 | p2) / (p3 | p4) / p5
  violindata <- sce_merge_res$sce@meta.data %>%
                                   select(.,orig.ident,nCount_RNA,nFeature_RNA,percent.mt) %>% 
                                   melt(.,id.vars="orig.ident") 
  p6 <- ggplot(data=violindata,mapping=aes(x=orig.ident,y=value,fill=orig.ident))+
         geom_violin()+
         facet_wrap(vars(variable),scales = "free_y") +
         theme_bw()+
         theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  (p1 | p2) / (p3 | p4) / p6   
}

#天津肿瘤P2_656098样本数据质控
{
  sce_P2_656098 <- FastSeurat(count=NULL,
                    dir.name="/Users/biofly/project/shijian/data/tianjian_ICC/P2_656098/mRNA/outs/filtered_feature_bc_matrix",
                    obj=NULL,#min.cells=3,
                    min.features=200,
                    max.features=7500,
                    percent.mt.num=50,
                    plot=F,
                    pcSelect=30,
                    project="P2_656098",
                    nfeatures=2000,
                    vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                    all.scale = F,
                    npcs = 50,
                    resolution = 0.8,
                    harmony=F,
                    doublet=T,
                    isMarkers=T,
                    cellCycle=T,
                    perplexity = 30,
                    features=NULL,
                    filepath=NULL,
                    outdir="Results",
                    names="love")
  save(sce_P2_656098,file="./Results/P2_656098_seurat.rda")
  sce_P2_656098$initial_sce #8188个细胞
  sce_P2_656098$sce #3203个细胞 2963个单细胞 240个双细胞
  print(sce_P2_656098$v) #min.features=200,max.features=2500,percent.mt.num=3
  print(sce_P2_656098$e) #pcSelect=30
  
  sce_P2_656098$sce #3203个细胞 2963个单细胞 240个双细胞
  print(sce_P2_656098$v) #min.features=200,max.features=4000,percent.mt.num=50
  print(sce_P2_656098$e) #pcSelect=30
  
  p1 <- sce_P2_656098$v
  p2 <- sce_P2_656098$e
  p3 <- DimPlot(sce_P2_656098$sce, reduction = "umap",label.size = 4,repel=T,label=T)
  p4 <- VlnPlot(sce_P2_656098$sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), idents=NULL,ncol = 3,pt.size=0)
  p2+p1+p3+p4+plot_layout(nrow=2)
  
  library(dplyr)
  #sce$sce.markers <- FindAllMarkers(object = sce$sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
  cl_markers <- sce_P2_656098$sce.markers
  top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  dh <- DoHeatmap(sce_P2_656098$sce, features = top10_cl_markers$gene) + NoLegend()
}

#天津肿瘤样本合并
{
  sce_merge <- FastMergeSeurat(objList=list(sce$sce,sce_P2_656098$sce),
                               object.names=c("P1_653303","P2_656098"),
                               project.name.default="ICC")
  sce_merge_res <- FastSeurat(count=NULL,
                              dir.name=NULL,
                              obj=sce_merge,#min.cells=3,
                              min.features=200,
                              max.features=7500,
                              percent.mt.num=50,
                              plot=F,
                              pcSelect=30,
                              project="ICC",
                              nfeatures=2000,
                              vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                              all.scale = F,
                              npcs = 50,
                              resolution = 0.5,
                              harmony=F,
                              doublet=F,
                              perplexity = 30,
                              isMarkers=T,
                              cellCycle=T,
                              features=NULL,
                              filepath=NULL,
                              outdir="Results",
                              names="love")
  saveRDS(sce_merge_res,file="./Results/P1_P2_sce_merge_noHamony_res.rds")
  sce_merge_res$sce@reductions$harmony@cell.embeddings
  sce_merge_res$sce@meta.data
  p1 <- sce_merge_res$v
  #p4 <- VlnPlot(sce_merge_res$sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), idents=NULL,ncol = 3,pt.size=0)
  head(sce_merge_res$sce@meta.data)
  bardata <- sce_merge_res$sce@meta.data %>% 
    group_by(orig.ident,Doublet) %>% 
    summarise(num=n())
  p1 <- ggplot(data=bardata,mapping=aes(x=orig.ident,y=num,fill=Doublet)) +
    geom_bar(stat = "identity") +
    #geom_text_repel() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  p2 <- sce_merge_res$e
  p3 <- DimPlot(sce_merge_res$sce, reduction = "umap",group.by="seurat_clusters",label.size = 4,repel=T,label=T)
  p4 <- DimPlot(sce_merge_res$sce, reduction = "umap",group.by="orig.ident",label.size = 4,repel=T,label=T)
  p5 <- VlnPlot(sce_merge_res$sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), idents=NULL,ncol = 3,pt.size=0)
  (p1 | p2) / (p3 | p4) / p5
  violindata <- sce_merge_res$sce@meta.data %>%
    select(.,orig.ident,nCount_RNA,nFeature_RNA,percent.mt) %>% 
    melt(.,id.vars="orig.ident") 
  p6 <- ggplot(data=violindata,mapping=aes(x=orig.ident,y=value,fill=orig.ident))+
    geom_violin()+
    facet_wrap(vars(variable),scales = "free_y") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  (p1 | p2) / (p3 | p4) / p6   
  
  library(dplyr)
  #sce$sce.markers <- FindAllMarkers(object = sce$sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
  cl_markers <- sce_merge_res$sce.markers
  top5_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  dh <- DoHeatmap(sce_merge_res$sce, features = top5_cl_markers$gene) + NoLegend()
  
  p2 <- VlnPlot(sce_merge_res$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p3 <- VlnPlot(sce_merge_res$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p4 <- VlnPlot(sce_merge_res$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  pdf("./Figure/compined plots.pdf",width = 8, height = 8)
  CombinePlots(plots = list(p2,p3,p4),nrow=3)
  dev.off()
}

#P3_651038
{
  sce_P3_651038 <- FastSeurat(count=NULL,
                              dir.name="/Volumes/HDD3/tianjin_ICC/P3_651038/mRNA/filtered_feature_bc_matrix",
                              obj=NULL,#min.cells=3,
                              min.features=0,
                              max.features=300000,
                              percent.mt.num=100,
                              species=c("human","mouse")[1],
                              plot=F,
                              pcSelect=30,
                              project="P3_651038",
                              nfeatures=2000,
                              vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                              all.scale = F,
                              npcs = 50,
                              resolution = 0.5,
                              harmony=F,
                              doublet=F,
                              isMarkers=F,
                              cellCycle=F,
                              perplexity = 30,
                              features=NULL,
                              filepath=NULL,
                              outdir="Results",
                              names="love")
  sce_P3_651038$initial_sce #8188个细胞
  sce_P3_651038$sce #3203个细胞 2963个单细胞 240个双细胞
  print(sce_P3_651038$v) #min.features=200,max.features=2500,percent.mt.num=3
  print(sce_P3_651038$e) #pcSelect=30
  sce_P3_651038$sce #3203个细胞 2963个单细胞 240个双细胞
  print(sce_P3_651038$v) #min.features=200,max.features=4000,percent.mt.num=50
  print(sce_P3_651038$e) #pcSelect=30
  p1 <- sce_P3_651038$v
  p2 <- sce_P3_651038$e
  p2 <- DimPlot(sce_P3_651038$sce, reduction = "umap",group.by="seurat_clusters",label.size = 6,repel=T,label=T)
  p3 <- DimPlot(sce_P3_651038$sce, reduction = "umap",group.by="orig.ident",label.size = 6,repel=T,label=T)
  p2 | p3
  p4 <- VlnPlot(sce_P3_651038$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce_P3_651038$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce_P3_651038$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  
  sce_P3_651038 <- FastSeurat(count=NULL,
                              dir.name="/Volumes/HDD3/tianjin_ICC/P3_651038/mRNA/filtered_feature_bc_matrix",
                              obj=NULL,#min.cells=3,
                              species=c("human","mouse")[1],
                              min.features=200,
                              max.features=10000,
                              percent.mt.num=50,
                              plot=F,
                              pcSelect=30,
                              project="P3_651038",
                              nfeatures=2000,
                              vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                              all.scale = F,
                              npcs = 50,
                              resolution = 0.5,
                              harmony=F,
                              doublet=T,
                              isMarkers=T,
                              cellCycle=T,
                              perplexity = 30,
                              features=NULL,
                              filepath=NULL,
                              outdir="Results",
                              names="love")
  saveRDS(sce_P3_651038,file="./Results/data/P3_651038_seurat.rds")
  p1 <- sce_P3_651038$v
  p1
  bardata <- sce_P3_651038$sce@meta.data %>% 
    group_by(orig.ident) %>% 
    summarise(num=n())
  p0 <- ggplot(data=bardata,mapping=aes(x=orig.ident,y=num,fill=orig.ident)) +
    geom_bar(stat = "identity") +
    #geom_text_repel() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  p2 <- DimPlot(sce_P3_651038$sce, reduction = "umap",group.by="seurat_clusters",label.size = 6,repel=T,label=T)
  p3 <- DimPlot(sce_P3_651038$sce, reduction = "umap",group.by="orig.ident",label.size = 6,repel=T,label=T)
  p2 | p3
  p4 <- VlnPlot(sce_P3_651038$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce_P3_651038$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce_P3_651038$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  
  cl_markers <- sce_P3_651038$sce.markers
  top20_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  top20_cl_markers %>% filter(.,cluster==9) -> cluster9_top20_cl_markers
  top20_cl_markers %>% filter(.,cluster==11) -> cluster11_top20_cl_markers
  dh <- DoHeatmap(sce_P3_651038$sce, features = cluster9_top20_cl_markers$gene,slot = "data") + NoLegend()
  dh <- DoHeatmap(sce_P3_651038$sce, features = cluster11_top20_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  
  dh <- DoHeatmap(sce_P3_651038$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
  top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  dh <- DoHeatmap(sce_P3_651038$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
}

#P4_662693
{
  sce_P4_662693 <- FastSeurat(count=NULL,
                              dir.name="/Volumes/HDD3/tianjin_ICC/P4_662693/mRNA/filtered_feature_bc_matrix",
                              obj=NULL,#min.cells=3,
                              min.features=0,
                              max.features=300000,
                              percent.mt.num=100,
                              species=c("human","mouse")[1],
                              plot=F,
                              pcSelect=30,
                              project="P4_662693",
                              nfeatures=2000,
                              vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                              all.scale = F,
                              npcs = 50,
                              resolution = 0.5,
                              harmony=F,
                              doublet=F,
                              isMarkers=F,
                              cellCycle=F,
                              perplexity = 30,
                              features=NULL,
                              filepath=NULL,
                              outdir="Results",
                              names="love")
  sce_P4_662693$initial_sce #8188个细胞
  sce_P3_651038$sce #3203个细胞 2963个单细胞 240个双细胞
  print(sce_P4_662693$v) #min.features=200,max.features=2500,percent.mt.num=3
  print(sce_P4_662693$e) #pcSelect=30
  sce_P4_662693$sce #3203个细胞 2963个单细胞 240个双细胞
  print(sce_P4_662693$v) #min.features=200,max.features=4000,percent.mt.num=50
  print(sce_P4_662693$e) #pcSelect=30
  p1 <- sce_P4_662693$v
  p2 <- sce_P4_662693$e
  p3 <- DimPlot(sce_P4_662693$sce, reduction = "umap",label.size = 4,repel=T,label=T)
  
  p2 <- DimPlot(sce_P4_662693$sce, reduction = "umap",group.by="seurat_clusters",label.size = 6,repel=T,label=T)
  p3 <- DimPlot(sce_P4_662693$sce, reduction = "umap",group.by="orig.ident",label.size = 6,repel=T,label=T)
  p2 | p3
  p4 <- VlnPlot(sce_P4_662693$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce_P4_662693$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce_P4_662693$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  
  sce_P4_662693 <- FastSeurat(count=NULL,
                              dir.name="/Volumes/HDD3/tianjin_ICC/P4_662693/mRNA/filtered_feature_bc_matrix",
                              obj=NULL,#min.cells=3,
                              species=c("human","mouse")[1],
                              min.features=100,
                              max.features=7500,
                              percent.mt.num=50,
                              plot=F,
                              pcSelect=30,
                              project="P4_662693",
                              nfeatures=2000,
                              vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                              all.scale = F,
                              npcs = 50,
                              resolution = 0.5,
                              harmony=F,
                              doublet=T,
                              isMarkers=T,
                              cellCycle=T,
                              perplexity = 30,
                              features=NULL,
                              filepath=NULL,
                              outdir="Results",
                              names="love")
  saveRDS(sce_P4_662693,file="./Results/data/P4_662693_seurat.rds")
  p1 <- sce_P4_662693$v
  p1
  p2 <- DimPlot(sce_P4_662693$sce, reduction = "umap",group.by="seurat_clusters",label.size = 6,repel=T,label=T)
  p2 
  p4 <- VlnPlot(sce_P4_662693$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce_P4_662693$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce_P4_662693$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  #marker热图
  cl_markers <- sce_P4_662693$sce.markers
  cl_markers %>%filter(.,cluster==4) %>% top_n(n = 20, wt = avg_log2FC) -> cluster4_top10_cl_markers
  cl_markers %>%filter(.,cluster==15) %>% top_n(n = 20, wt = avg_log2FC) -> cluster15_top10_cl_markers
  cl_markers %>%filter(.,cluster==11) %>% top_n(n = 20, wt = avg_log2FC) -> cluster11_top10_cl_markers
  dh <- DoHeatmap(sce_P4_662693$sce, features = cluster11_top10_cl_markers$gene,slot = "data") + NoLegend()
  dh <- DoHeatmap(sce_P4_662693$sce, features = cluster15_top10_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  print(dh)
  
  top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  dh <- DoHeatmap(sce_P4_662693$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
}

#P1 P2 P3 P4 merge 
{
  P1_653303_seurat <- get(load(file="./Results/data/P1_653303_seurat.rda"))
  P2_656098_seurat <- get(load(file="./Results/data/P2_656098_seurat.rda"))
  P3_651038_seurat <- readRDS(file="./Results/data/P3_651038_seurat.rds")
  P4_662693_seurat <- readRDS(file="./Results/data/P4_662693_seurat.rds")
  sce_merge <- FastMergeSeurat(objList=list(P1_653303_seurat$sce,P2_656098_seurat$sce,P3_651038_seurat$sce,P4_662693_seurat$sce),
                               object.names=c("P1_653303","P2_656098","P3_651038","P4_662693"),
                               species = c("human","mouse")[1],
                               project.name.default="ICC")
  sce_merge_res <- FastSeurat(count=NULL,
                              dir.name=NULL,
                              obj=sce_merge,#min.cells=3,
                              species=c("human","mouse")[1],
                              min.features=0,
                              max.features=300000,
                              percent.mt.num=100,
                              plot=F,
                              pcSelect=30,
                              project="ICC",
                              nfeatures=2000,
                              vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                              all.scale = F,
                              npcs = 50,
                              resolution = 0.5,
                              harmony=F,
                              doublet=F,
                              perplexity = 30,
                              isMarkers=F,
                              cellCycle=F,
                              features=NULL,
                              filepath=NULL,
                              outdir="Results",
                              names="love")
  saveRDS(sce_merge_res,file="./Results/data/sce_merge_res.rds")
  bardata <- sce_merge_res$sce@meta.data %>% 
    dplyr::group_by(orig.ident) %>% 
    dplyr::summarise(num=n())
  p0 <- ggplot(data=bardata,mapping=aes(x=orig.ident,y=num,fill=orig.ident)) +
    geom_bar(stat = "identity") +
    #geom_text_repel() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  bardata <- sce_merge_res$sce@meta.data %>% 
    dplyr::group_by(orig.ident) %>% 
    dplyr::summarise(num=n())
  VlnPlot(sce_merge_res$sce,features = c('nCount_RNA','nFeature_RNA','percent.mt'),group.by = "orig.ident",pt.size=0)
  p1 <- sce_merge_res$v
  p1
  p2 <- DimPlot(sce_merge_res$sce, reduction = "umap",group.by="seurat_clusters",label.size = 6,repel=T,label=T)
  p3 <- DimPlot(sce_merge_res$sce, reduction = "umap",group.by="orig.ident",label.size = 6,repel=T,label=T)
  p2 | p3 
  p4 <- VlnPlot(sce_merge_res$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce_merge_res$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce_merge_res$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  #marker热图
  sce_merge_res <- readRDS("./Results/data/sce_merge_res.rds")
  cl_markers <- sce_merge_res$sce.markers
  top10_cl_markers <- cl_markers %>% dplyr::group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  dh <- DoHeatmap(sce_merge_res$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  ################----------------------------####################
  #因为观察到亚群之间的marker有很多线粒体基因和免疫球蛋白基因
  #所以对合并后的4个样本去除其线粒体基因，核糖体基因后再重新聚类
  sce_4_merge <- FastSeurat_1(count=NULL,
                              dir.name=NULL,
                              obj=sce_merge,
                              species=c("human","mouse")[1],
                              min.features=0,
                              max.features=300000,
                              percent.mt.num=100,
                              plot=F,
                              pcSelect=30,
                              project="ICC",
                              nfeatures=2000,
                              vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                              all.scale = F,
                              npcs = 50,
                              resolution = 0.5,
                              harmony=F,
                              doublet=F,
                              isMarkers=F,
                              rmOtherGene=T,
                              perplexity = 30,
                              cellCycle=T,
                              features=NULL,
                              filepath=NULL,
                              outdir="Results",
                              names="love")
  
  saveRDS(sce_4_merge,file = "./Results/data/sce_4_merge.rds")
  bardata <- sce_4_merge$sce@meta.data %>% 
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
  VlnPlot(sce_4_merge$sce,features = c('nCount_RNA','nFeature_RNA','percent.mt',"percent.ribo"),group.by = "orig.ident",pt.size=0)
  head(sce_4_merge$sce@meta.data)
  p1 <- dittoDimPlot(sce_4_merge$sce, reduction.use  = "umap",var="orig.ident",size = 1,do.label = TRUE, do.ellipse = TRUE,legend.size = 9,shape.legend.size=9,labels.size=8,do.raster=T,raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 18))
  p2 <- dittoDimPlot(sce_4_merge$sce, reduction.use  = "umap",var="seurat_clusters",size = 1,do.label = TRUE, do.ellipse = TRUE,legend.size = 9,shape.legend.size=9,labels.size=8,do.raster=T,raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 18))
  p3 <- plot_grid(p1,p2)
  ggsave(p3,file="")
  p4 <- VlnPlot(sce_4_merge$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce_4_merge$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce_4_merge$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  
  FeaturePlot_scCustom(sce_4_merge$sce, features = EPCAM, order = F)
  FeaturePlot_scCustom(sce_4_merge$sce, features = cholangiocytes, order = F)
  FeaturePlot_scCustom(sce_4_merge$sce, features = hepatocytes,order = F)
  FeaturePlot_scCustom(sce_4_merge$sce, features = fibroblasts,order = F)
  FeaturePlot_scCustom(sce_4_merge$sce, features = endothelial,order = F)
  FeaturePlot_scCustom(sce_4_merge$sce, features = All.Immune,order = F)
  FeaturePlot_scCustom(sce_4_merge$sce, features = B.cell,order = F)
  FeaturePlot_scCustom(sce_4_merge$sce, features = plasma,order = F)
  FeaturePlot_scCustom(sce_4_merge$sce, features = T.cell,order = F)
  FeaturePlot_scCustom(sce_4_merge$sce, features = proliferative,order = F)
  FeaturePlot_scCustom(sce_4_merge$sce, features = NK,order = F)
  FeaturePlot_scCustom(sce_4_merge$sce, features = Dendritic.cell,order = F)
  FeaturePlot_scCustom(sce_4_merge$sce, features = Macrophage,order = F)
  FeaturePlot_scCustom(sce_4_merge$sce, features = Mast,order = F)
  FeaturePlot_scCustom(sce_4_merge$sce, features = Monocyte,order = F)
  FeaturePlot_scCustom(sce_4_merge$sce, features = pDC,order = F)
  #观察不确定亚群19,25,6的marker基因
  cl_markers <- FindAllMarkers(sce_4_merge$sce) 
  top10_cl_markers <- cl_markers %>% dplyr::group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
  dh <- DoHeatmap(sce_4_merge$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  
  NK_T_sce <- subset(sce_4_merge$sce, seurat_clusters %in% c(0,1,4,12,21))
  NK_T_sce <- FastSeurat(count=NULL,
                            dir.name=NULL,
                            obj=NK_T_sce,
                            species=c("human","mouse")[1],
                            min.features=0,
                            max.features=300000,
                            percent.mt.num=100,
                            plot=F,
                            pcSelect=30,
                            project="ICC",
                            nfeatures=2000,
                            vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                            all.scale = F,
                            npcs = 50,
                            resolution = 0.5,
                            harmony=F,
                            doublet=F,
                            isMarkers=T,
                            rmOtherGene=T,
                            perplexity = 30,
                            cellCycle=T,
                            features=NULL,
                            filepath=NULL,
                            outdir="Results",
                            names="love")
  saveRDS(NK_T_sce,file = "./Results/data/NK_T_sce.rds")
  bardata <- NK_T_sce$sce@meta.data %>% 
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
  VlnPlot(NK_T_sce$sce,features = c('nCount_RNA','nFeature_RNA','percent.mt',"percent.ribo"),group.by = "orig.ident",pt.size=0)
  p1 <- dittoDimPlot(NK_T_sce$sce, reduction.use  = "umap",var="orig.ident",size = 1,do.label = TRUE, do.ellipse = TRUE,legend.size = 9,shape.legend.size=9,labels.size=8,do.raster=T,raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 18))
  p2 <- dittoDimPlot(NK_T_sce$sce, reduction.use  = "umap",var="seurat_clusters",size = 1,do.label = TRUE, do.ellipse = TRUE,legend.size = 9,shape.legend.size=9,labels.size=8,do.raster=T,raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 18))
  p3 <- plot_grid(p1,p2)
  p4 <- VlnPlot(NK_T_sce$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(NK_T_sce$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(NK_T_sce$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  Treg <- c("FOXP3")
  Exhausted_T <- c("LAG3")
  Proliferative_T <- c("TOP2A")
  FeaturePlot_scCustom(NK_T_sce$sce, features = c(T.cell,Treg),order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = EPCAM, order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = cholangiocytes, order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = hepatocytes,order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = fibroblasts,order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = endothelial,order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = All.Immune,order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = B.cell,order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = plasma,order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = T.cell,order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = Treg,order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = Exhausted_T,order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = Proliferative_T,order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = NK,order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = Dendritic.cell,order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = Macrophage,order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = Mast,order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = Monocyte,order = F)
  FeaturePlot_scCustom(NK_T_sce$sce, features = pDC,order = F)
  

  cl_markers <- NK_T_sce$sce.markers
  top10_cl_markers <- cl_markers %>% dplyr::group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
  dh <- DoHeatmap(NK_T_sce$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  
  Myeloid_sce <- subset(sce_4_merge$sce, seurat_clusters %in% c(7,11))
  Myeloid_sce <- FastSeurat(count=NULL,
                         dir.name=NULL,
                         obj=Myeloid_sce,
                         species=c("human","mouse")[1],
                         min.features=0,
                         max.features=300000,
                         percent.mt.num=100,
                         plot=F,
                         pcSelect=30,
                         project="ICC",
                         nfeatures=2000,
                         vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                         all.scale = F,
                         npcs = 50,
                         resolution = 0.5,
                         harmony=F,
                         doublet=F,
                         isMarkers=T,
                         rmOtherGene=T,
                         perplexity = 30,
                         cellCycle=T,
                         features=NULL,
                         filepath=NULL,
                         outdir="Results",
                         names="love") 
  save(Myeloid_sce,file = "./Results/data/Myeloid_sce.rds")
  bardata <- Myeloid_sce$sce@meta.data %>% 
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
  VlnPlot(Myeloid_sce$sce,features = c('nCount_RNA','nFeature_RNA','percent.mt',"percent.ribo"),group.by = "orig.ident",pt.size=0)
  p1 <- dittoDimPlot(Myeloid_sce$sce, reduction.use  = "umap",var="orig.ident",size = 1,do.label = TRUE, do.ellipse = TRUE,legend.size = 9,shape.legend.size=9,labels.size=8,do.raster=T,raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 18))
  p2 <- dittoDimPlot(Myeloid_sce$sce, reduction.use  = "umap",var="seurat_clusters",size = 1,do.label = TRUE, do.ellipse = TRUE,legend.size = 9,shape.legend.size=9,labels.size=8,do.raster=T,raster.dpi = 500) + theme(legend.text = element_text(face="bold",size = 18))
  p3 <- plot_grid(p1,p2)
  p4 <- VlnPlot(Myeloid_sce$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(Myeloid_sce$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(Myeloid_sce$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  
  FeaturePlot_scCustom(Myeloid_sce$sce, features = c(T.cell,Treg),order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = EPCAM, order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = cholangiocytes, order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = hepatocytes,order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = fibroblasts,order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = endothelial,order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = All.Immune,order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = B.cell,order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = plasma,order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = T.cell,order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = Treg,order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = Exhausted_T,order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = Proliferative_T,order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = NK,order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = Dendritic.cell,order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = Macrophage,order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = Mast,order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = Monocyte,order = F)
  FeaturePlot_scCustom(Myeloid_sce$sce, features = pDC,order = F)
  
  
  
} 

#国人ICC
{
  library(data.table)
  sampleAnno <- fread("/Users/biofly/project/shijian/data/GSE138709/filereport_read_run_PRJNA576876_tsv.txt")
  GSM4116579_ajacent <- fread("/Users/biofly/project/shijian/data/GSE138709/GSE138709_RAW/GSM4116579_ICC_18_Adjacent_UMI.csv.gz")
  GSM4116586_ajacent <- fread("/Users/biofly/project/shijian/data/GSE138709/GSE138709_RAW/GSM4116586_ICC_25_Adjacent_UMI.csv.gz")
  GSM4116582_ajacent <- fread("/Users/biofly/project/shijian/data/GSE138709/GSE138709_RAW/GSM4116582_ICC_23_Adjacent_UMI.csv.gz")
  GSM4116580 <- fread("/Users/biofly/project/shijian/data/GSE138709/GSE138709_RAW/GSM4116580_ICC_18_Tumor_UMI.csv.gz")
  GSM4116581 <- fread("/Users/biofly/project/shijian/data/GSE138709/GSE138709_RAW/GSM4116581_ICC_20_Tumor_UMI.csv.gz")
  GSM4116583 <- fread("/Users/biofly/project/shijian/data/GSE138709/GSE138709_RAW/GSM4116583_ICC_23_Tumor_UMI.csv.gz")
  GSM4116584 <- fread("/Users/biofly/project/shijian/data/GSE138709/GSE138709_RAW/GSM4116584_ICC_24_Tumor1_UMI.csv.gz")
  GSM4116585 <- fread("/Users/biofly/project/shijian/data/GSE138709/GSE138709_RAW/GSM4116585_ICC_24_Tumor2_UMI.csv.gz")
  GSM4116579_ajacent <- as.data.frame(GSM4116579_ajacent)
  GSM4116586_ajacent <- as.data.frame(GSM4116586_ajacent)
  GSM4116582_ajacent <- as.data.frame(GSM4116582_ajacent)
  GSM4116580 <- as.data.frame(GSM4116580)
  GSM4116581 <- as.data.frame(GSM4116581)
  GSM4116583 <- as.data.frame(GSM4116583)
  GSM4116584 <- as.data.frame(GSM4116584)
  GSM4116585 <- as.data.frame(GSM4116585)
  rownames(GSM4116579_ajacent) <- GSM4116579_ajacent$V1
  rownames(GSM4116586_ajacent) <- GSM4116586_ajacent$V1
  rownames(GSM4116582_ajacent) <- GSM4116582_ajacent$V1
  rownames(GSM4116580) <- GSM4116580$V1
  rownames(GSM4116581) <- GSM4116581$V1
  rownames(GSM4116583) <- GSM4116583$V1
  rownames(GSM4116584) <- GSM4116584$V1
  rownames(GSM4116585) <- GSM4116585$V1
  GSM4116579_ajacent <- GSM4116579_ajacent[,-1]
  GSM4116586_ajacent <- GSM4116586_ajacent[,-1]
  GSM4116582_ajacent <- GSM4116582_ajacent[,-1]
  GSM4116580 <- GSM4116580[,-1]
  GSM4116581 <- GSM4116581[,-1]
  GSM4116583 <- GSM4116583[,-1]
  GSM4116584 <- GSM4116584[,-1]
  GSM4116585 <- GSM4116585[,-1]
  GSM4116579_ajacent <- as.matrix(GSM4116579_ajacent)  
  GSM4116586_ajacent <- as.matrix(GSM4116586_ajacent)
  GSM4116582_ajacent <- as.matrix(GSM4116582_ajacent)
  GSM4116580 <- as.matrix(GSM4116580)
  GSM4116581 <- as.matrix(GSM4116581)
  GSM4116583 <- as.matrix(GSM4116583)
  GSM4116584 <- as.matrix(GSM4116584)
  GSM4116585 <- as.matrix(GSM4116585)
  expr.list <- list(GSM4116579=GSM4116579_ajacent,GSM4116586=GSM4116586_ajacent,GSM4116582=GSM4116582_ajacent,GSM4116580=GSM4116580,GSM4116581=GSM4116581,
                    GSM4116583=GSM4116583,GSM4116584=GSM4116584,GSM4116585=GSM4116585)
  sceList <- list()
  for(i in 1:length(expr.list)){
    print(i)
    tt <- FastSeurat(count=expr.list[[i]],
                     dir.name=NULL,
                     obj=NULL,#min.cells=3,
                     species=c("human","mouse")[1],
                     min.features=0,
                     max.features=300000,
                     percent.mt.num=50,
                     plot=F,
                     pcSelect=30,
                     project=names(expr.list)[i],
                     nfeatures=2000,
                     vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                     all.scale = F,
                     npcs = 50,
                     resolution = 0.5,
                     harmony=F,
                     doublet=F,
                     perplexity = 30,
                     isMarkers=F,
                     cellCycle=F,
                     features=NULL,
                     filepath=NULL,
                     outdir="Results",
                     names="love")  
    sceList[[i]] <- tt
  }
  {
    sceList[[7]]$v
    p2 <- DimPlot(sceList[[7]]$sce, reduction = "umap",group.by="seurat_clusters",label.size = 6,repel=T,label=T)
    p3 <- DimPlot(sceList[[7]]$sce, reduction = "umap",group.by="orig.ident",label.size = 6,repel=T,label=T)
    p2 | p3
    p4 <- VlnPlot(sceList[[7]]$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
      geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
    p5 <- VlnPlot(sceList[[7]]$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
      geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
    p6 <- VlnPlot(sceList[[7]]$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
      geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
      theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
    CombinePlots(plots = list(p4,p5,p6),nrow=3)
  }
  names(sceList) <- names(expr.list)
  sce.list <- lapply(sceList,function(x){
    x$sce
  })
  P1_653303_seurat <- get(load(file="./Results/data/P1_653303_seurat.rda"))
  P2_656098_seurat <- get(load(file="./Results/data/P2_656098_seurat.rda"))
  P3_651038_seurat <- readRDS(file="./Results/data/P3_651038_seurat.rds")
  P4_662693_seurat <- readRDS(file="./Results/data/P4_662693_seurat.rds")
  # sce_merge_res <- readRDS(file="./Results/data/sce_merge_res.rds")
  sce_merge <- FastMergeSeurat(objList=c(P1_653303_seurat$sce,P2_656098_seurat$sce,P3_651038_seurat$sce,P4_662693_seurat$sce,sce.list),
                               object.names=c("P1_653303","P2_656098","P3_651038","P4_662693",names(expr.list)),
                               species = c("human","mouse")[1],
                               project.name.default="ICC")
  unique(sce_merge@meta.data$orig.ident)
  sce_merge <- FastSeurat(count=NULL,
                     dir.name=NULL,
                     obj=sce_merge,#min.cells=3,
                     species=c("human","mouse")[1],
                     min.features=0,
                     max.features=300000,
                     percent.mt.num=50,
                     plot=F,
                     pcSelect=30,
                     project=names(expr.list)[i],
                     nfeatures=2000,
                     vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                     all.scale = F,
                     npcs = 50,
                     resolution = 0.5,
                     harmony=T,
                     doublet=F,
                     perplexity = 30,
                     isMarkers=T,
                     cellCycle=T,
                     features=NULL,
                     filepath=NULL,
                     outdir="Results",
                     names="love") 
  saveRDS(sce_merge,file="./Results/data/sce_merge_paper.rds")
  bardata <- sce_merge$sce@meta.data %>% 
    dplyr::group_by(orig.ident) %>% 
    dplyr::summarise(num=n())
  p0 <- ggplot(data=bardata,mapping=aes(x=orig.ident,y=num,fill=orig.ident)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  VlnPlot(sce_merge$sce,features = c('nCount_RNA','nFeature_RNA','percent.mt'),group.by = "orig.ident",pt.size=0)
  p1 <- sce_merge$v
  p1
  p2 <- DimPlot(sce_merge$sce, reduction = "umap",group.by="seurat_clusters",label.size = 6,repel=T,label=T,cols = zylcolor40)
  p3 <- DimPlot(sce_merge$sce, reduction = "umap",group.by="orig.ident",label.size = 6,repel=T,label=T,cols = zylcolor40)
  p2 | p3 
  p4 <- VlnPlot(sce_merge$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce_merge$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce_merge$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  #marker热图
  cl_markers <- sce_merge$sce.markers
  top10_cl_markers <- cl_markers %>% dplyr::group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  dh <- DoHeatmap(sce_merge$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
}

#我们的4个样本精细分群以及精细去除双细胞
{
  sce_merge_res  <- readRDS(file="./Results/data/sce_merge_res.rds")
  immune.cell <- c(21,12,6,14,27,22,0,2,4,7,13) 
  T.cell <- c(0,2,7,13,21)
  NK <- 4
  B.cell <- c(14,27,22)
  Myeloid <- c(6,12)
  epcam <- c(26,24,8,16,15,10,23,19,3,18,9,1,5,11)
  head(sce_merge_res$sce@meta.data)
  sce_T_NK <- subset(sce_merge_res$sce,seurat_clusters %in% c(0,2,7,13,21,4))
  sce_B <- subset(sce_merge_res$sce,seurat_clusters %in% c(14,27,22))
  sce_Myeloid <- subset(sce_merge_res$sce,seurat_clusters %in% c(6,12))
  sce_EPCAM <- subset(sce_merge_res$sce,seurat_clusters %in% c(26,24,8,16,15,10,23,19,3,18,9,1,5,11))
  sce_B_res <- FastSeurat(count=NULL,
                              dir.name=NULL,
                              obj=sce_B,#min.cells=3,
                              species=c("human","mouse")[1],
                              min.features=0,
                              max.features=300000,
                              percent.mt.num=100,
                              plot=F,
                              pcSelect=30,
                              project="ICC",
                              nfeatures=2000,
                              vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                              all.scale = F,
                              npcs = 50,
                              resolution = 0.5,
                              harmony=F,
                              doublet=F,
                              perplexity = 30,
                              isMarkers=T,
                              cellCycle=T,
                              features=NULL,
                              filepath=NULL,
                              outdir="Results",
                              names="love")
  saveRDS(sce_B_res,file="./Results/data/sce_B_res.rds")       
  head(sce_B_res$sce@meta.data)
  p2 <- DimPlot(sce_B_res$sce, reduction = "umap",group.by="seurat_clusters",label.size = 6,repel=T,label=T)
  p3 <- DimPlot(sce_B_res$sce, reduction = "umap",group.by="orig.ident",label.size = 6,repel=T,label=T)
  p2 | p3 
  p4 <- VlnPlot(sce_B_res$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce_B_res$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce_B_res$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  # marker着色
  aDC <- "LAMP3"
  EPCAM <- c("EPCAM","KRT19","KRT7") #上皮细胞
  cholangiocytes <- c("FYXD2","TM4SF4", "ANXA4") #胆管细胞
  hepatocytes  <- c("APOC3","FABP1","APOA1") # 肝脏细胞
  fibroblasts <- c("ACTA2" ,"COL1A2") #成纤维细胞
  endothelial <- c("ENG","VWF") #内皮细胞
  All.Immune <- "PTPRC" #免疫细胞
  B.cell <- c("BANK1","CD19","CD79A","CD79B","MS4A1")
  T.cell <- c("CD3D","CD3E","CD8A","CD4")
  CD8.T <- c("CD8A","GZMK","GZMB","PIK3R1","TUBA4A","TNFRSF9","TOP2A","MKI67")
  CD4.T <- c("CD4","IL7R")
  NK.cell <- c("KLRF1")
  Treg.cell <- c("FOXP3")
  Exhausted.cell <- c("LAG3","TIGIT","TIM3","HAVCR1")
  Dendritic.cell <- c("CLEC9A","CD1C","THBD","ITGAX")
  Macrophage <- c("CD68","CD14","HLA-DRA","HLA-DRB1","APOE","MMP9","CTSK")
  Mast <- c("KIT","TPSAB1","TPSB2","ENPP3")
  Monocyte <- c("LYZ","CD68","CD163","ITGAX","S100A8","S100A9")
  NK <- c("FCGR3A","FCGR3B","NCAM1","PVRIG","TIGIT","CD7","FGFBP2","KLRF1")
  pDC <- c("LILRA4","UGCG","CLEC4C")
  plasma <- c("SDC1","CD38","SLAMF7","DERL3","MZB1","IGHA1","IGHA2")
  proliferative <- c("MKI67","PCNA","TK1","TOP2A","TUBB","TYMS") #增殖
  FeaturePlot(sce_B_res$sce, features = NK.cell)
  FeaturePlot(sce_B_res$sce, features = Treg.cell)
  FeaturePlot(sce_B_res$sce, features = Exhausted.cell)
  FeaturePlot(sce_B_res$sce, features = CD4.T)
  FeaturePlot(sce_B_res$sce, features = aDC)
  FeaturePlot(sce_B_res$sce, features = EPCAM)
  FeaturePlot(sce_B_res$sce, features = cholangiocytes)
  FeaturePlot(sce_B_res$sce, features = hepatocytes)
  FeaturePlot(sce_B_res$sce, features = fibroblasts)
  FeaturePlot(sce_B_res$sce, features = endothelial)
  FeaturePlot(sce_B_res$sce, features = All.Immune)
  FeaturePlot(sce_B_res$sce, features = B.cell)
  FeaturePlot(sce_B_res$sce, features = plasma)
  FeaturePlot(sce_B_res$sce, features = T.cell)
  FeaturePlot(sce_B_res$sce, features = NK)
  FeaturePlot(sce_B_res$sce, features = Dendritic.cell)
  FeaturePlot(sce_B_res$sce, features = Macrophage)
  FeaturePlot(sce_B_res$sce, features = Mast)
  FeaturePlot(sce_B_res$sce, features = Monocyte)
  FeaturePlot(sce_B_res$sce, features = pDC)
  FeaturePlot(sce_B_res$sce, features = proliferative)
  cl_markers <- sce_B_res$sce.markers
  top10_cl_markers <- cl_markers %>% dplyr::group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  dh <- DoHeatmap(sce_B_res$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  B.doublet <- c(21,13,9)
  head(sce_B_res$sce@meta.data)
  sce_B_single_res <- subset(sce_B_res$sce,seurat_clusters != 21 & seurat_clusters!= 13 & seurat_clusters != 9)
  head(sce_B_res$sce@meta.data)
  sce_EPCAM_res <- FastSeurat(count=NULL,
                                dir.name=NULL,
                                obj=sce_EPCAM,#min.cells=3,
                                species=c("human","mouse")[1],
                                min.features=0,
                                max.features=300000,
                                percent.mt.num=100,
                                plot=F,
                                pcSelect=30,
                                project="ICC",
                                nfeatures=2000,
                                vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                                all.scale = F,
                                npcs = 50,
                                resolution = 0.5,
                                harmony=F,
                                doublet=F,
                                perplexity = 30,
                                isMarkers=T,
                                cellCycle=T,
                                features=NULL,
                                filepath=NULL,
                                outdir="Results",
                                names="love")
  saveRDS(sce_EPCAM_res,file="./Results/data/sce_EPCAM_res.rds")       
  head(sce_EPCAM_res$sce@meta.data)
  p2 <- DimPlot(sce_EPCAM_res$sce, reduction = "umap",group.by="seurat_clusters",label.size = 6,repel=T,label=T)
  p3 <- DimPlot(sce_EPCAM_res$sce, reduction = "umap",group.by="orig.ident",label.size = 6,repel=T,label=T)
  p2 | p3 
  p4 <- VlnPlot(sce_EPCAM_res$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce_EPCAM_res$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce_EPCAM_res$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  # marker着色
  aDC <- "LAMP3"
  EPCAM <- c("EPCAM","KRT19","KRT7") #上皮细胞
  cholangiocytes <- c("FYXD2","TM4SF4", "ANXA4") #胆管细胞
  hepatocytes  <- c("APOC3","FABP1","APOA1") # 肝脏细胞
  fibroblasts <- c("ACTA2" ,"COL1A2") #成纤维细胞
  endothelial <- c("ENG","VWF") #内皮细胞
  All.Immune <- "PTPRC" #免疫细胞
  B.cell <- c("BANK1","CD19","CD79A","CD79B","MS4A1")
  T.cell <- c("CD3D","CD3E","CD8A","CD4")
  CD8.T <- c("CD8A","GZMK","GZMB","PIK3R1","TUBA4A","TNFRSF9","TOP2A","MKI67")
  CD4.T <- c("CD4","IL7R")
  NK.cell <- c("KLRF1")
  Treg.cell <- c("FOXP3")
  Exhausted.cell <- c("LAG3","TIGIT","TIM3","HAVCR1")
  Dendritic.cell <- c("CLEC9A","CD1C","THBD","ITGAX")
  Macrophage <- c("CD68","CD14","HLA-DRA","HLA-DRB1","APOE","MMP9","CTSK")
  Mast <- c("KIT","TPSAB1","TPSB2","ENPP3")
  Monocyte <- c("LYZ","CD68","CD163","ITGAX","S100A8","S100A9")
  NK <- c("FCGR3A","FCGR3B","NCAM1","PVRIG","TIGIT","CD7","FGFBP2","KLRF1")
  pDC <- c("LILRA4","UGCG","CLEC4C")
  plasma <- c("SDC1","CD38","SLAMF7","DERL3","MZB1","IGHA1","IGHA2")
  proliferative <- c("MKI67","PCNA","TK1","TOP2A","TUBB","TYMS") #增殖
  FeaturePlot(sce_EPCAM_res$sce, features = NK.cell)
  FeaturePlot(sce_EPCAM_res$sce, features = Treg.cell)
  FeaturePlot(sce_EPCAM_res$sce, features = Exhausted.cell)
  FeaturePlot(sce_EPCAM_res$sce, features = CD4.T)
  FeaturePlot(sce_EPCAM_res$sce, features = aDC)
  FeaturePlot(sce_EPCAM_res$sce, features = EPCAM)
  FeaturePlot(sce_EPCAM_res$sce, features = cholangiocytes)
  FeaturePlot(sce_EPCAM_res$sce, features = hepatocytes)
  FeaturePlot(sce_EPCAM_res$sce, features = fibroblasts)
  FeaturePlot(sce_EPCAM_res$sce, features = endothelial)
  FeaturePlot(sce_EPCAM_res$sce, features = All.Immune)
  FeaturePlot(sce_EPCAM_res$sce, features = B.cell)
  FeaturePlot(sce_EPCAM_res$sce, features = plasma)
  FeaturePlot(sce_EPCAM_res$sce, features = T.cell)
  FeaturePlot(sce_EPCAM_res$sce, features = NK)
  FeaturePlot(sce_EPCAM_res$sce, features = Dendritic.cell)
  FeaturePlot(sce_EPCAM_res$sce, features = Macrophage)
  FeaturePlot(sce_EPCAM_res$sce, features = Mast)
  FeaturePlot(sce_EPCAM_res$sce, features = Monocyte)
  FeaturePlot(sce_EPCAM_res$sce, features = pDC)
  FeaturePlot(sce_EPCAM_res$sce, features = proliferative)
  cl_markers <- sce_EPCAM_res$sce.markers
  top10_cl_markers <- cl_markers %>% dplyr::group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  dh <- DoHeatmap(sce_EPCAM_res$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  epcam.doublet <- c(21,13,9)
  head(sce_EPCAM_res$sce@meta.data)
  sce_EPCAM_single_res <- subset(sce_EPCAM_res$sce,seurat_clusters != 21 & seurat_clusters!= 13 & seurat_clusters != 9)
  head(sce_Myeloid_res$sce@meta.data)
  sce_Myeloid_res <- FastSeurat(count=NULL,
                             dir.name=NULL,
                             obj=sce_Myeloid,#min.cells=3,
                             species=c("human","mouse")[1],
                             min.features=0,
                             max.features=300000,
                             percent.mt.num=100,
                             plot=F,
                             pcSelect=30,
                             project="ICC",
                             nfeatures=2000,
                             vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                             all.scale = F,
                             npcs = 50,
                             resolution = 0.5,
                             harmony=F,
                             doublet=F,
                             perplexity = 30,
                             isMarkers=T,
                             cellCycle=T,
                             features=NULL,
                             filepath=NULL,
                             outdir="Results",
                             names="love")
  saveRDS(sce_Myeloid_res,file="./Results/data/sce_Myeloid_res.rds")
  head(sce_Myeloid_res$sce@meta.data)
  p2 <- DimPlot(sce_Myeloid_res$sce, reduction = "umap",group.by="seurat_clusters",label.size = 6,repel=T,label=T)
  p3 <- DimPlot(sce_Myeloid_res$sce, reduction = "umap",group.by="orig.ident",label.size = 6,repel=T,label=T)
  p2 | p3 
  p4 <- VlnPlot(sce_Myeloid_res$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce_Myeloid_res$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce_Myeloid_res$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  # marker着色
  aDC <- "LAMP3"
  EPCAM <- c("EPCAM","KRT19","KRT7") #上皮细胞
  cholangiocytes <- c("FYXD2","TM4SF4", "ANXA4") #胆管细胞
  hepatocytes  <- c("APOC3","FABP1","APOA1") # 肝脏细胞
  fibroblasts <- c("ACTA2" ,"COL1A2") #成纤维细胞
  endothelial <- c("ENG","VWF") #内皮细胞
  All.Immune <- "PTPRC" #免疫细胞
  B.cell <- c("BANK1","CD19","CD79A","CD79B","MS4A1")
  T.cell <- c("CD3D","CD3E","CD8A","CD4")
  CD8.T <- c("CD8A","GZMK","GZMB","PIK3R1","TUBA4A","TNFRSF9","TOP2A","MKI67")
  CD4.T <- c("CD4","IL7R")
  NK.cell <- c("KLRF1")
  Treg.cell <- c("FOXP3")
  Exhausted.cell <- c("LAG3","TIGIT","TIM3","HAVCR1")
  Dendritic.cell <- c("CLEC9A","CD1C","THBD","ITGAX")
  Macrophage <- c("CD68","CD14","HLA-DRA","HLA-DRB1","APOE","MMP9","CTSK")
  Mast <- c("KIT","TPSAB1","TPSB2","ENPP3")
  Monocyte <- c("LYZ","CD68","CD163","ITGAX","S100A8","S100A9")
  NK <- c("FCGR3A","FCGR3B","NCAM1","PVRIG","TIGIT","CD7","FGFBP2","KLRF1")
  pDC <- c("LILRA4","UGCG","CLEC4C")
  plasma <- c("SDC1","CD38","SLAMF7","DERL3","MZB1","IGHA1","IGHA2")
  proliferative <- c("MKI67","PCNA","TK1","TOP2A","TUBB","TYMS") #增殖
  FeaturePlot(sce_Myeloid_res$sce, features = NK.cell)
  FeaturePlot(sce_Myeloid_res$sce, features = Treg.cell)
  FeaturePlot(sce_Myeloid_res$sce, features = Exhausted.cell)
  FeaturePlot(sce_Myeloid_res$sce, features = CD4.T)
  FeaturePlot(sce_Myeloid_res$sce, features = aDC)
  FeaturePlot(sce_Myeloid_res$sce, features = EPCAM)
  FeaturePlot(sce_Myeloid_res$sce, features = cholangiocytes)
  FeaturePlot(sce_Myeloid_res$sce, features = hepatocytes)
  FeaturePlot(sce_Myeloid_res$sce, features = fibroblasts)
  FeaturePlot(sce_Myeloid_res$sce, features = endothelial)
  FeaturePlot(sce_Myeloid_res$sce, features = All.Immune)
  FeaturePlot(sce_Myeloid_res$sce, features = B.cell)
  FeaturePlot(sce_Myeloid_res$sce, features = plasma)
  FeaturePlot(sce_Myeloid_res$sce, features = T.cell)
  FeaturePlot(sce_Myeloid_res$sce, features = NK)
  FeaturePlot(sce_Myeloid_res$sce, features = Dendritic.cell)
  FeaturePlot(sce_Myeloid_res$sce, features = Macrophage)
  FeaturePlot(sce_Myeloid_res$sce, features = Mast)
  FeaturePlot(sce_Myeloid_res$sce, features = Monocyte)
  FeaturePlot(sce_Myeloid_res$sce, features = pDC)
  FeaturePlot(sce_Myeloid_res$sce, features = proliferative)
  cl_markers <- sce_Myeloid_res$sce.markers
  top10_cl_markers <- cl_markers %>% dplyr::group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  dh <- DoHeatmap(sce_Myeloid_res$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
  sce_T_NK_res <- FastSeurat(count=NULL,
                             dir.name=NULL,
                             obj=sce_T_NK,#min.cells=3,
                             species=c("human","mouse")[1],
                             min.features=0,
                             max.features=300000,
                             percent.mt.num=100,
                             plot=F,
                             pcSelect=30,
                             project="ICC",
                             nfeatures=2000,
                             vars.to.regress = NULL,  #c("nCounts_RNA", "percent.mito")
                             all.scale = F,
                             npcs = 50,
                             resolution = 0.5,
                             harmony=F,
                             doublet=F,
                             perplexity = 30,
                             isMarkers=T,
                             cellCycle=T,
                             features=NULL,
                             filepath=NULL,
                             outdir="Results",
                             names="love")
  saveRDS(sce_T_NK_res,file="./Results/data/sce_T_NK_res.rds")
  sce_T_NK_res <- readRDS(file="./Results/data/sce_T_NK_res.rds")
  head(sce_T_NK_res$sce@meta.data)
  p2 <- DimPlot(sce_T_NK_res$sce, reduction = "umap",group.by="seurat_clusters",label.size = 6,repel=T,label=T)
  p3 <- DimPlot(sce_T_NK_res$sce, reduction = "umap",group.by="orig.ident",label.size = 6,repel=T,label=T)
  p2 | p3 
  p4 <- VlnPlot(sce_T_NK_res$sce,features='nCount_RNA',ncol=1,pt.size=0,log=T,same.y.lims=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p5 <- VlnPlot(sce_T_NK_res$sce,features='nFeature_RNA',ncol=1,pt.size=0,log=T) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")
  p6 <- VlnPlot(sce_T_NK_res$sce,features='percent.mt',ncol=1,pt.size=0,y.max = 100) +
    geom_boxplot(outlier.size=0,width=0.3,show.legend=F,notchwidth = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 1),legend.position = "none")
  CombinePlots(plots = list(p4,p5,p6),nrow=3)
  # marker着色
  aDC <- "LAMP3"
  EPCAM <- c("EPCAM","KRT19","KRT7") #上皮细胞
  cholangiocytes <- c("FYXD2","TM4SF4", "ANXA4") #胆管细胞
  hepatocytes  <- c("APOC3","FABP1","APOA1") # 肝脏细胞
  fibroblasts <- c("ACTA2" ,"COL1A2") #成纤维细胞
  endothelial <- c("ENG","VWF") #内皮细胞
  All.Immune <- "PTPRC" #免疫细胞
  B.cell <- c("BANK1","CD19","CD79A","CD79B","MS4A1")
  T.cell <- c("CD3D","CD3E","CD8A","CD4")
  CD8.T <- c("CD8A","GZMK","GZMB","PIK3R1","TUBA4A","TNFRSF9","TOP2A","MKI67")
  CD4.T <- c("CD4","IL7R")
  NK.cell <- c("KLRF1")
  Treg.cell <- c("FOXP3")
  Exhausted.cell <- c("LAG3","TIGIT","TIM3","HAVCR1")
  Dendritic.cell <- c("CLEC9A","CD1C","THBD","ITGAX")
  Macrophage <- c("CD68","CD14","HLA-DRA","HLA-DRB1","APOE","MMP9","CTSK")
  Mast <- c("KIT","TPSAB1","TPSB2","ENPP3")
  Monocyte <- c("LYZ","CD68","CD163","ITGAX","S100A8","S100A9")
  NK <- c("FCGR3A","FCGR3B","NCAM1","PVRIG","TIGIT","CD7","FGFBP2","KLRF1")
  pDC <- c("LILRA4","UGCG","CLEC4C")
  plasma <- c("SDC1","CD38","SLAMF7","DERL3","MZB1","IGHA1","IGHA2")
  proliferative <- c("MKI67","PCNA","TK1","TOP2A","TUBB","TYMS") #增殖
  FeaturePlot(sce_T_NK_res$sce, features = NK.cell)
  FeaturePlot(sce_T_NK_res$sce, features = Treg.cell)
  FeaturePlot(sce_T_NK_res$sce, features = Exhausted.cell)
  FeaturePlot(sce_T_NK_res$sce, features = CD4.T)
  FeaturePlot(sce_T_NK_res$sce, features = aDC)
  FeaturePlot(sce_T_NK_res$sce, features = EPCAM)
  FeaturePlot(sce_T_NK_res$sce, features = cholangiocytes)
  FeaturePlot(sce_T_NK_res$sce, features = hepatocytes)
  FeaturePlot(sce_T_NK_res$sce, features = fibroblasts)
  FeaturePlot(sce_T_NK_res$sce, features = endothelial)
  FeaturePlot(sce_T_NK_res$sce, features = All.Immune)
  FeaturePlot(sce_T_NK_res$sce, features = B.cell)
  FeaturePlot(sce_T_NK_res$sce, features = plasma)
  FeaturePlot(sce_T_NK_res$sce, features = T.cell)
  FeaturePlot(sce_T_NK_res$sce, features = NK)
  FeaturePlot(sce_T_NK_res$sce, features = Dendritic.cell)
  FeaturePlot(sce_T_NK_res$sce, features = Macrophage)
  FeaturePlot(sce_T_NK_res$sce, features = Mast)
  FeaturePlot(sce_T_NK_res$sce, features = Monocyte)
  FeaturePlot(sce_T_NK_res$sce, features = pDC)
  FeaturePlot(sce_T_NK_res$sce, features = proliferative)
  cl_markers <- sce_T_NK_res$sce.markers
  top10_cl_markers <- cl_markers %>% dplyr::group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  dh <- DoHeatmap(sce_T_NK_res$sce, features = top10_cl_markers$gene,slot = "data") + NoLegend()
  print(dh)
}

