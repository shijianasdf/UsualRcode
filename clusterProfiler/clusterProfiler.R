library(clusterProfiler)
library(org.Hs.eg.db)
x <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2", 
       "ORM1",  "IGFBP1", "PTHLH", "GPC3",  "IGFBP3","TOB1",    "MITF",   "NDRG1", 
       "NR1H4", "FGFR3",  "PVR",   "IL6",   "PTPRM", "ERBB2",   "NID2",   "LAMB1", 
       "COMP",  "PLS3",   "MCAM",  "SPP1",  "LAMC1", "COL4A2",  "COL4A1", "MYOC",  
       "ANXA4", "TFPI2",  "CST6",  "SLPI",  "TIMP2", "CPM",     "GGT1",   "NNMT",
       "MAL",   "EEF1A2", "HGD",   "TCN2",  "CDA",   "PCCA",    "CRYM",   "PDXK",  
       "STC1",  "WARS",  "HMOX1", "FXYD2", "RBP4",   "SLC6A12", "KDELR3", "ITM2B");
#id转换
geneList<- bitr(x,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL"),org.Hs.eg.db);
#GO富集分析
go <- enrichGO(geneList$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,qvalueCutoff = 0.2,keyType = 'ENTREZID');
# 进行简单的可视化
barplot(go,showCategory=20,drop=T)
dotplot(go,showCategory=50)
#KEGG富集分析
kegg <- enrichKEGG(geneList$ENTREZID, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE);
# 进行简单的可视化
dotplot(kegg, showCategory=30);

#GSEA
coad_diff_expression <- get(load("D:/R/Project/Learning/Transcription/COAD_filter_DESeq_normalized_results.RData"))
pos <- order(coad_diff_expression$DESeq_res$log2FoldChange,decreasing =T)
temp <- cbind.data.frame(rownames(coad_diff_expression$DESeq_res)[pos],coad_diff_expression$DESeq_res$log2FoldChange[pos]);
colnames(temp) <- c("SYMBOL","LOG2FOLDCHANGE")
temp$SYMBOL <- as.character(temp$SYMBOL)
SYMBOL_ENTREZID <- select(org.Hs.eg.db,keys=temp$SYMBOL,keytype="SYMBOL",columns="ENTREZID")
merge_table <- merge(temp,SYMBOL_ENTREZID,by="SYMBOL")
pos1 <- order(merge_table$LOG2FOLDCHANGE,decreasing =T)
merge_table <- merge_table[pos1,]
merge_table <- na.omit(merge_table)
geneList <- merge_table$LOG2FOLDCHANGE
names(geneList) <- as.character(merge_table$ENTREZID)
#mget(temp$SYMBOL, org.Hs.egSYMBOL2EG, ifnotfound=NA)
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)
gsea <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE, pvalueCutoff = 0.8);

#GSEA KEGG分析
gseKEGG <- gseKEGG(geneList, organism = "hsa", keyType = "kegg", exponent = 1,
  nPerm = 1000, minGSSize = 10, maxGSSize = 500,
  pvalueCutoff = 0.5, pAdjustMethod = "BH", verbose = TRUE,
  use_internal_data = FALSE, seed = FALSE, by = "fgsea")
gseaplot(gseKEGG,geneSetID=rownames(gseKEGG[1,])) 

library(pathview)
pathview(gene.data = geneList, pathway.id = 'hsa04658',species="hsa", limit=list(gene=max(abs(geneList)), cpd=1))