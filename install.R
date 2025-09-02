#设置国内镜像
options()$repos 
options()$BioC_mirror
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options()$repos 
options()$BioC_mirror

#如果下载不下来，就去网站下载下来，本地安装
#' 包的安装
install.packages("devtools","remotes","BiocManager")
install.packages(c('tibble','plyr','pacman','ggplot2','ggpubr','RColorBrewer','tuneR','readr','progress','dplyr','data.table','magrittr','ggplotify',
                   'stringr','reshape2','tidyr','ggsignif','ggrepel','ggthemes','ggsci','ggExtra','ggforce','cowplot','patchwork'))
install.packages('WGCNA','e1071', 'preprocessCore',"ppcor",  "timeROC", "pracma",'survival', 'survminer',"caret",'pROC','umap')
install.packages(c("FactoMineR", "factoextra","XML","glmnet","maxstat","rvest","caret","pheatmap","ComplexHeatmap","circlize","ggdendro","ggtree","igraph",
                   "ggraph","tidygraph","tidytext","tidyverse","plotly","shiny","shinyjs","shinyBS","shinyWidgets",
                   "shinycssloaders","shinythemes","DT","visNetwork","networkD3","DiagrammeR",
                   "rmarkdown","knitr"))
install.packages(c("ggtreeExtra","ggtree","ggtreeAssist","ggtreeExtra","ggtreeify","ggtreeLayout","ggtreeify","ggtreeDendro","ggtreeGrob"))
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
devtools::install_github("IOBR/IOBR", build_vignettes=TRUE)
devtools::install_github("GreenleafLab/chromVARmotifs")
devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
devtools::install_github("immunogenomics/presto")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
devtools::install_github("immunogenomics/harmony", build_vignettes=TRUE)
#下载到服务器，然后用install_local()来安装
devtools::install_local("PATH/TO/DIRECTORY/CytoTRACE_0.3.3.tar.gz")

# https://bioconductor.org/packages/release/bioc/html/GEOquery.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("KEGG.db","dittoSeq","scCustomize"),ask = F,update = F)
BiocManager::install(c("GSEABase","GSVA","clusterProfiler","GenomicRanges","GenomicFeatures","GenomicAlignments"),ask = F,update = F)
BiocManager::install(c("GEOquery","limma","impute" ),ask = F,update = F)
BiocManager::install(c("genefu","org.Hs.eg.db","hgu133plus2.db" ),ask = F,update = F)
BiocManager::install(c('airway','DESeq2','edgeR'),ask = F,update = F)
BiocManager::install(c('scran'),ask = F,update = F)
BiocManager::install(c("TxDb.Mmusculus.UCSC.mm10.knownGene","org.Mm.eg.db","TxDb.Hsapiens.UCSC.hg38.knownGene"),ask = F,update = F)
BiocManager::install(c('Seurat','monocle'),ask = F,update = F)
BiocManager::install(c('destiny','scRNAseq','dbscan','M3Drop','flexclust','mcclust'),ask = F,update = F)
BiocManager::install(c("biomaRt","sva","GO.db"),ask = F,update = F)
BiocManager::install(c('monocle3','Signac',"cicero","JASPAR2020","TFBSTools",
                       "BSgenome.Hsapiens.UCSC.hg38","EnsDb.Hsapiens.v86","chromVAR",
                       "ChIPseeker","SingleCellExperiment","NMF","scater","SC3"),ask = F,update = F)
BiocManager::install(c("motifmatchr","harmony"),ask = F,update = F)

library(ArchR)
library(Signac)
library(cicero)
library(data.table)
library(dplyr)
library(Matrix)
library(ggplot2)
library(scCustomize)
library(Seurat)
library(patchwork)
library(dittoSeq)
library(scales)
library(umap)
library(factoextra)
library(FactoMineR)
library(RColorBrewer)
library(JASPAR2020)
library(TFBSTools)
library(chromVARmotifs)
data("human_pwms_v2")
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(chromVAR)
library(rstatix)
library(org.Hs.eg.db)
library(data.table)
library(dittoSeq)
library(scCustomize)
library(Seurat)
library(CytoTRACE)
library(CytoTRACE2)
library(stringr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(tidytext)
library(GenomicRanges)
library(GenomicFeatures)
library(foreach)
library(doParallel)
library(data.table)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(SC3)
library(future)
library(patchwork)
library(cowplot)
library(ggsci)
library(NMF)
library(SingleCellExperiment)
library(scater)






