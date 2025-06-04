#' 包的安装
install.packages("devtools")
install.packages(c('plyr','pacman','ggplot2','ggpubr','RColorBrewer','tuneR','readr','progress','dplyr'))
install.packages('WGCNA')
install.packages(c("FactoMineR", "factoextra","XML","glmnet","maxstat","rvest"))
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)

options()$repos 
options()$BioC_mirror
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options()$repos 
options()$BioC_mirror
# https://bioconductor.org/packages/release/bioc/html/GEOquery.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("KEGG.db",ask = F,update = F)
BiocManager::install(c("GSEABase","GSVA","clusterProfiler" ),ask = F,update = F)
BiocManager::install(c("GEOquery","limma","impute" ),ask = F,update = F)
BiocManager::install(c("genefu","org.Hs.eg.db","hgu133plus2.db" ),ask = F,update = F)
BiocManager::install(c('airway','DESeq2','edgeR'),ask = F,update = F)
BiocManager::install(c('scran'),ask = F,update = F)
BiocManager::install(c("TxDb.Mmusculus.UCSC.mm10.knownGene","org.Mm.eg.db","TxDb.Hsapiens.UCSC.hg38.knownGene"),ask = F,update = F)
BiocManager::install(c('Seurat','monocle'),ask = F,update = F)
BiocManager::install(c('destiny','scRNAseq','dbscan','M3Drop','flexclust','mcclust'),ask = F,update = F)
BiocManager::install(c("biomaRt","sva","GO.db"),ask = F,update = F)
