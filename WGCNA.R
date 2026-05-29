#下载芯片表达数据
FastDownloadGEO(GSEid="GSE5281",filepath="D:",origin=F)
#加载表达数据，表型数据以及探针注释信息
GSE5281_expr <- fread("D:/GSE5281/GSE5281-GPL570-expression.MATRIX.txt")
GSE5281_pData <- fread("D:/GSE5281/GSE5281-GPL570-pData1.txt")
GSE5281_ProbeAnnotation <- fread("D:/GSE5281/GSE5281-GPL570-ProbeAnnotation.txt")
GSE5281_expr <- GSE5281_expr %>% tibble::column_to_rownames("V1")
GSE5281_pData <- as.data.frame(GSE5281_pData)
rownames(GSE5281_pData) <- GSE5281_pData$geo_accession
colnames(GSE5281_ProbeAnnotation)[11] <- "Gene.Symbol"
GSE5281_ProbeAnnotation <- as.data.frame(GSE5281_ProbeAnnotation)
load("E:/rProj/sCcRCC/data/common.annot.genecode.V48.2025.11.16.rda")
probeAnnotation_res <- probeAnnotation(exprSet=GSE5281_expr,probeAnnotation=GSE5281_ProbeAnnotation,id="ID",symbol="Gene.Symbol",
                                       entrezid="ENTREZ_GENE_ID",index="Gene.Symbol",
                                       fromtype=c("ENTREZID","SYMBOL","ENSEMBL")[2],
                                       totype = c("entrez","symbol","ensemble")[2])
probeAnnotation_res$exprSet[1:4,1:4]
probeAnnotation_res$exprSet <- log2eset(probeAnnotation_res$exprSet)
#res <- BeforeWGCNA(expr.matrix=as.matrix(probeAnnotation_res$exprSet),design=as.data.frame(GSE5281_pData),log.convert=F,mad.portion=0.75,mad.min = 0.01)
GSE5281_pData$Age <- stringr::str_split(GSE5281_pData$Age," ",simplify = T)[,1]
GSE5281_pData$Age[which(GSE5281_pData$Age == ">90")] <- "90"
GSE5281_pData$Age <- as.numeric(GSE5281_pData$Age)
GSE5281_pData$Disease_State <- ifelse(grepl("normal",GSE5281_pData$Disease_State),"normal","Alzheimer")
GSE5281_pData$Sex <- ifelse(grepl("female",GSE5281_pData$Sex),"female","male")
# design <- GetWGCNATrait(design=GSE5281_pData[,c("Disease_State","Sex","Age")],
#                         convert.list=list(Disease_State=list("Alzheimer" = 1,"normal" = 0),
#                                           Sex=list("male"= 1,"female"= 0)))
limma_res <- limma1(counts=probeAnnotation_res$exprSet,design=GSE5281_pData,contrast.col="Disease_State",method = c("voom","trend")[2],
                   contrast.level = c("Alzheimer","normal"), # NP相对N0计算差异,NO是control组
                   cutoff.lFC = 0.5,cutoff.padj = 0.05,save.file = T,
                   report  = T,db=common.annot,names = "Alzheimer")
load("D:/GSE5281/Alzheimer_limmaList.rda")
test_res <- compareGenes_Boxplot(expr=probeAnnotation_res$exprSet,design=GSE5281_pData,
                     features=c("ATP5C1","PSMB3","ACACB","SLC35E1"),
                                 sampleCol="geo_accession",conditionCol="Disease_State")
test_res$p  
intersect(limma_res$siggenes,c("ACP2","ADK","CDK10","DHFR","MAPKAPK3","MGST3","PKN3","THNSL2"))  

#筛选的四个基因与已知阿尔默兹海默症关键基因的相关性
cor_res <- MakeCorrelationMatrix(x=t(probeAnnotation_res$exprSet), y=t(probeAnnotation_res$exprSet), size= 500, 
                                 adjust= "BH",use=c("everything", "all.obs", "complete.obs", "na.or.complete","pairwise.complete.obs")[1],
                                 method = c("pearson", "kendall", "spearman")[1])
save(cor_res,file="D:/GSE5281/cor_genes.rda")
cor_res <- get(load(file="D:/GSE5281/cor_genes.rda"))
cor_res1 <- outputDF(lst=cor_res,symmetric= TRUE)
ob_genes <- c("ACP2","ADK","CDK10","DHFR")
Alzheimer_genes <- c("APOE","TREM2", "ABCA7","APP", "PSEN1", "PSEN2", "BIN1", "SORL1","CLU", "CR1", "CD33", "MS4A4A", "MS4A6A", "TREML4", "TREM3", 
                     "INPP5D", "ABCG1", "ABCA1", "CASS4", "CELF1", "FERMT2", "HLA-DRB5", "HLA-DRB1", "NME8", "PTK2B", "SLC24A4", "ZCWPW1", "MEF2C")
head(cor_res1)
filter_cor <- cor_res1 %>% filter(.,((V1 %in% ob_genes) & (V2 %in% Alzheimer_genes)) | ((V1 %in% Alzheimer_genes) & (V2 %in% ob_genes))) %>% 
                           filter(.,r > 0.4 & p < 0.05)
#INPP5D    DHFR  0.349223988 5.614436e-06 7.179241e-05

plotCor.data <- t(probeAnnotation_res$exprSet[c("ADK","CDK10","ACP2","DHFR","ABCA1","PTK2B","PSEN1","ABCA7","MEF2C","SORL1"),])
FastGgpairs(data=plotCor.data, columns=c("ADK","CDK10","ACP2","DHFR","ABCA1","PTK2B","PSEN1","ABCA7","MEF2C","SORL1"),title="test")

#识别的基因和阿尔海默症基因相关性
plotCor.data %<>% as.data.frame() %>% tibble::rownames_to_column("sample") %>% 
                 reshape2::melt(.,id.vars=c("sample","ACP2","ADK","CDK10","DHFR"),variable.name= "Alzheimer_genes",value.name = "expression") %>% 
                 reshape2::melt(.,id.vars=c("sample","Alzheimer_genes","expression"),variable.name= "ob_genes",value.name = "ob_expression")
p <- ggplot(data=as.data.frame(plotCor.data),mapping=aes(x=ob_expression,y=expression)) +
  geom_point(alpha=.3,size=3) + 
  geom_smooth(method = "lm", formula = y ~ x, se = T) + 
  stat_cor(method = "pearson", label.x.npc = 'left', label.y.npc = 'top') +
  facet_grid(Alzheimer_genes~ob_genes,scales = "free_y",labeller = label_wrap_gen(width=10)) + 
  #geom_text_repel(mapping=aes(label=!!sym(sample.col))) +
  theme_bw() +
  theme(text = element_text(size=20),strip.text = element_text(size = 10,face = "bold"),
        panel.margin = unit(0, "lines"),strip.placement = "outside")

temp.plotCor.data <- c()
for(i in 1:nrow(plotCor.data)){
  for(j in 1:nrow(filter_cor)){
    if( ((plotCor.data[i,2] == filter_cor[j,1]) & (plotCor.data[i,4] == filter_cor[j,2])) | ((plotCor.data[i,4] == filter_cor[j,1]) & (plotCor.data[i,2] == filter_cor[j,2]))){
      temp.plotCor.data <- rbind(temp.plotCor.data,
                                 data.frame(sample=plotCor.data$sample[i],
                                            Alzheimer_genes=plotCor.data$Alzheimer_genes[i],
                                            expression=plotCor.data$expression[i],
                                            ob_genes=plotCor.data$ob_genes[i],
                                            ob_expression=plotCor.data$ob_expression[i]))
    }
  }
}
temp.plotCor.data %<>% mutate(Alzheimer_genes = as.character(Alzheimer_genes),
                             ob_genes = as.character(ob_genes))
p <- ggplot(data=as.data.frame(temp.plotCor.data),mapping=aes(x=ob_expression,y=expression)) +
  geom_point(alpha=.3,size=3) + 
  geom_smooth(method = "lm", formula = y ~ x, se = T) + 
  stat_cor(method = "pearson", label.x.npc = 'left', label.y.npc = 'top') +
  facet_grid(Alzheimer_genes~ob_genes,scales = "free_y",labeller = label_wrap_gen(width=10)) + 
  #geom_text_repel(mapping=aes(label=!!sym(sample.col))) +
  theme_bw() +
  theme(text = element_text(size=20),strip.text = element_text(size = 10,face = "bold"),
        panel.margin = unit(0, "lines"),strip.placement = "outside") 

plotCor.data1 <- t(probeAnnotation_res$exprSet[c("ADK","CDK10","ACP2","DHFR","ABCA1","PTK2B","PSEN1","ABCA7","MEF2C","SORL1"),])
plotCor.data1 %<>% as.data.frame() %>% tibble::rownames_to_column("sample") 
ob_genes <- c("ACP2","ADK","CDK10","DHFR")
p.list <- c()
for(i in 1:nrow(filter_cor)){
  if(filter_cor$V1[i] %in% ob_genes){
    p <- FastCorScatterPlot(plotCor.data1,x=as.character(filter_cor$V1[i]),y=as.character(filter_cor$V2[i]),facet_by=NULL,sample.col="sample",
                       method=c("pearson","spearman","kendall")[1],
                       label.x.npc=c('right', 'left', 'center', 'centre', 'middle')[2],
                       label.y.npc=c('bottom', 'top', 'center', 'centre', 'middle')[2])
  }else{
    p <- FastCorScatterPlot(plotCor.data1,x=as.character(filter_cor$V2[i]),y=as.character(filter_cor$V1[i]),facet_by=NULL,sample.col="sample",
                            method=c("pearson","spearman","kendall")[1],
                            label.x.npc=c('right', 'left', 'center', 'centre', 'middle')[2],
                            label.y.npc=c('bottom', 'top', 'center', 'centre', 'middle')[2])
  }
  p.list <- c(p.list,p)
}
p <- (p.list[[1]] + p.list[[2]] + p.list[[3]]) /  (p.list[[4]] + p.list[[5]] + p.list[[6]])
ggsave(p,filename="D:/GSE5281/Alzheimer_correation.pdf",width = 12,height=8)
  
WGCNA_res <- FastWGCNA(expr.matrix=probeAnnotation_res$exprSet,design=GSE5281_pData,log.convert = F,
                       check = T, mad.portion = 1, mad.min = 0,contrast.col = "Disease_State",contrast.control = "normal",
                       verbose = c(goodSamplesGenes = 3,
                        pickSoftThreshold = 5,
                        blockwiseModules = 3), maxBlockSize = 50000, minModuleSize = 30,
                       type = c("unsigned","signed","signed hybrid")[1],
                      cutoff.pval = 0.05,corType = c("pearson","bicor")[1],hub_cutoffSigGM=0.2,
                      hub_MM = 0.8,hub_WeightedQ = 0.01,parallel = T,report = F,two.step = F,
                      save.path = "WGCNA2",names = "XieJing")
save(probeAnnotation_res,GSE5281_pData,WGCNA_res,file="D:/GSE5281/WGCNA_res.rda")
load(file="D:/GSE5281/WGCNA_res.rda")
intersect(WGCNA_res$Hub$Block1$HubGenes,c("ACP2","ADK","CDK10","DHFR","MAPKAPK3","MGST3","PKN3","THNSL2")) 
intersect(colnames(WGCNA_res$dataExpr),c("ACP2","ADK","CDK10","DHFR","MAPKAPK3","MGST3","PKN3","THNSL2")) 
dim(WGCNA_res$dataExpr)
dim(probeAnnotation_res$exprSet)
WGCNA_res$sft$powerEstimate
WGCNA_res$net$MEs[1:4,1:4]
moduleColors <- labels2colors(WGCNA_res$net$colors)
module_df <- data.frame(cluster=unname(WGCNA_res$net$colors),gene=names(WGCNA_res$net$colors),moduleColors=moduleColors)
ob_genes <- c("ACP2","ADK","CDK10","DHFR")
cyt <- get(load("E:/rProj/sCcRCC/XieJing/WGCNA2/XieJing_Block1_Cytoscape Network.rda"))
rm(TOM);gc()
module_df[match(ob_genes,module_df$gene),]
#cluster  gene moduleColors
#4792       3  ACP2        brown
#5605       1   ADK    turquoise
#5203       7 CDK10        black
#4660       5  DHFR        green
#提取brown ACP2
brown_genes <- cyt$nodeData$nodeName[which(cyt$nodeData[,3] == "brown")]
head(cyt$edgeData)
brown_edgeData <- cyt$edgeData[(cyt$edgeData$fromNode %in% brown_genes) & (cyt$edgeData$toNode %in% brown_genes),]
brown_edgeData %<>% filter(.,fromNode == "ACP2" | toNode == "ACP2") %>% 
                   filter(.,weight > 0.05) %>% 
                   bind_rows(.,brown_edgeData[brown_edgeData$fromNode == "ACP2" & brown_edgeData$toNode == "PTK2B",])
brown_edgeData %>% filter(.,((fromNode == "ACP2") & (toNode %in% c("PTK2B"))) | ((toNode == "ACP2") & (fromNode %in% c("PTK2B"))))
unique(c(brown_edgeData$fromNode,brown_edgeData$toNode))
brown_nodes <- data.frame(name=unique(c(brown_edgeData$fromNode,brown_edgeData$toNode)),
                          type=rep("genes",length(unique(c(brown_edgeData$fromNode,brown_edgeData$toNode)))))
brown_nodes$type[which(brown_nodes$name == "ACP2")] = "risk_gene"
brown_network <- FastnetworkVisualization(edges=brown_edgeData[,1:3],nodes=brown_nodes,
                         layout=c("stress")[1],directed = F)
pdf("D:/GSE5281/brown_network.pdf",width=10,height=10)
print(brown_network$p)
dev.off()
ACP2_brown <- brown_nodes$name
ACP2_brown <- na.omit(convert(ACP2_brown,fromtype = "SYMBOL",totype = "ENTREZID"))
geneList <- t(WGCNA_res$dataExpr)[,1]
geneList_names <- convert(names(geneList),fromtype = "SYMBOL",totype = "ENTREZID")
pos <- which(!is.na(geneList_names))
geneList <- geneList[pos]
names(geneList) <- geneList_names[pos]
ACP2_brown_GOres <- FastGO(ACP2_brown,geneList,organism = c("human","mouse")[1],
                         default.universe = F,classlevel = 2:2,
                         OrgDb  = NULL,  #c("org.Mm.eg.db","org.Hs.eg.db")
                         keyType = NULL,pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,qvalueCutoff  = 0.05,
                         cnet.showCategory = 5,verbose = TRUE,save.path = "ACP2_brown_GO",names = "XieJing")


#提取turquoise
turquoise_genes <- cyt$nodeData$nodeName[which(cyt$nodeData[,3] == "turquoise")]
turquoise_edgeData <- cyt$edgeData[(cyt$edgeData$fromNode %in% turquoise_genes) & (cyt$edgeData$toNode %in% turquoise_genes),]
turquoise_edgeData1 <- turquoise_edgeData %>% filter(.,fromNode == "ADK" | toNode == "ADK") %>% 
                                             filter(.,weight > 0.13) 
turquoise_edgeData2 <- turquoise_edgeData %>% filter(.,((fromNode == "ADK") & (toNode %in% c("PSEN1","MEF2C","SORL1"))) | ((toNode == "ADK") & (fromNode %in% c("PSEN1","MEF2C","SORL1"))))
turquoise_edgeData3 <- turquoise_edgeData1 %>% dplyr::bind_rows(turquoise_edgeData2)
unique(c(turquoise_edgeData3$fromNode,turquoise_edgeData3$toNode))
turquoise_nodes <- data.frame(name=unique(c(turquoise_edgeData3$fromNode,turquoise_edgeData3$toNode)),
                          type=rep("genes",length(unique(c(turquoise_edgeData3$fromNode,turquoise_edgeData3$toNode)))))
turquoise_nodes$type[which(turquoise_nodes$name == "ADK")] = "risk_gene"
turquoise_network <- FastnetworkVisualization(edges=turquoise_edgeData3[,1:3],nodes=turquoise_nodes,
                                          layout=c("stress")[1],directed = F)
pdf("D:/GSE5281/turquoise_network.pdf",width=10,height=10)
print(turquoise_network$p)
dev.off()
ADK_turquoise <- turquoise_nodes$name
ADK_turquoise <- na.omit(convert(ADK_turquoise,fromtype = "SYMBOL",totype = "ENTREZID"))
geneList <- t(WGCNA_res$dataExpr)[,1]
geneList_names <- convert(names(geneList),fromtype = "SYMBOL",totype = "ENTREZID")
pos <- which(!is.na(geneList_names))
geneList <- geneList[pos]
names(geneList) <- geneList_names[pos]
ADK_turquoise_GOres <- FastGO(ADK_turquoise,geneList,organism = c("human","mouse")[1],
                           default.universe = F,classlevel = 2:2,
                           OrgDb  = NULL,  #c("org.Mm.eg.db","org.Hs.eg.db")
                           keyType = NULL,pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,qvalueCutoff  = 0.05,
                           cnet.showCategory = 5,verbose = TRUE,save.path = "ADK_turquoise_GO",names = "XieJing")

#提取black CDK10
black_genes <- cyt$nodeData$nodeName[which(cyt$nodeData[,3] == "black")]
black_edgeData <- cyt$edgeData[(cyt$edgeData$fromNode %in% black_genes) & (cyt$edgeData$toNode %in% black_genes),]
black_edgeData1 <- black_edgeData %>% filter(.,fromNode == "CDK10" | toNode == "CDK10") %>% 
                                      filter(.,weight > 0.02) 
black_edgeData2 <- black_edgeData %>% filter(.,((fromNode == "CDK10") & (toNode == "ABCA7")) | ((toNode == "CDK10") & (fromNode == "ABCA7")))
black_edgeData3 <- black_edgeData1 %>% dplyr::bind_rows(black_edgeData2)
unique(c(black_edgeData3$fromNode,black_edgeData3$toNode))
black_nodes <- data.frame(name=unique(c(black_edgeData3$fromNode,black_edgeData3$toNode)),
                              type=rep("genes",length(unique(c(black_edgeData3$fromNode,black_edgeData3$toNode)))))
black_nodes$type[which(black_nodes$name == "CDK10")] = "risk_gene"
black_network <- FastnetworkVisualization(edges=black_edgeData3[,1:3],nodes=black_nodes,
                                              layout=c("stress")[1],directed = F)
pdf("D:/GSE5281/black_network.pdf",width=10,height=10)
print(black_network$p)
dev.off()
CDK10_black <- black_nodes$name
CDK10_black <- na.omit(convert(CDK10_black,fromtype = "SYMBOL",totype = "ENTREZID"))
geneList <- t(WGCNA_res$dataExpr)[,1]
geneList_names <- convert(names(geneList),fromtype = "SYMBOL",totype = "ENTREZID")
pos <- which(!is.na(geneList_names))
geneList <- geneList[pos]
names(geneList) <- geneList_names[pos]
CDK10_black_GOres <- FastGO(CDK10_black,geneList,organism = c("human","mouse")[1],
                              default.universe = F,classlevel = 2:2,
                              OrgDb  = NULL,  #c("org.Mm.eg.db","org.Hs.eg.db")
                              keyType = NULL,pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,qvalueCutoff  = 0.05,
                              cnet.showCategory = 5,verbose = TRUE,save.path = "CDK10_black_GO",names = "XieJing")

#提取green DHFR
green_genes <- cyt$nodeData$nodeName[which(cyt$nodeData[,3] == "green")]
green_edgeData <- cyt$edgeData[(cyt$edgeData$fromNode %in% green_genes) & (cyt$edgeData$toNode %in% green_genes),]
green_edgeData1 <- green_edgeData %>% filter(.,fromNode == "DHFR" | toNode == "DHFR") %>% 
                                      filter(.,weight > 0.025) 
green_edgeData2 <- green_edgeData %>% filter(.,((fromNode == "DHFR") & (toNode == "ABCA1")) | ((toNode == "CDK10") & (fromNode == "ABCA1")))
green_edgeData3 <- green_edgeData1 %>% dplyr::bind_rows(green_edgeData2)
green_edgeData3 <- unique(green_edgeData3)
unique(c(green_edgeData3$fromNode,green_edgeData3$toNode))
green_nodes <- data.frame(name=unique(c(green_edgeData3$fromNode,green_edgeData3$toNode)),
                          type=rep("genes",length(unique(c(green_edgeData3$fromNode,green_edgeData3$toNode)))))
green_nodes$type[which(green_nodes$name == "DHFR")] = "risk_gene"
green_network <- FastnetworkVisualization(edges=green_edgeData3[,1:3],nodes=green_nodes,
                                          layout=c("stress")[1],directed = F)
pdf("D:/GSE5281/green_network.pdf",width=10,height=10)
print(green_network$p)
dev.off()
DHFR_green <- green_nodes$name
DHFR_green <- na.omit(convert(DHFR_green,fromtype = "SYMBOL",totype = "ENTREZID"))
geneList <- t(WGCNA_res$dataExpr)[,1]
geneList_names <- convert(names(geneList),fromtype = "SYMBOL",totype = "ENTREZID")
pos <- which(!is.na(geneList_names))
geneList <- geneList[pos]
names(geneList) <- geneList_names[pos]
DHFR_green_GOres <- FastGO(DHFR_green,geneList,organism = c("human","mouse")[1],
                            default.universe = F,classlevel = 2:2,
                            OrgDb  = NULL,  #c("org.Mm.eg.db","org.Hs.eg.db")
                            keyType = NULL,pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,qvalueCutoff  = 0.05,
                            cnet.showCategory = 5,verbose = TRUE,save.path = "DHFR_green_GO",names = "XieJing")



{
  head(WGCNA_res$net$colors,10)
  head(labels2colors(WGCNA_res$net$colors),10)
  ME_genes <- cbind(as.data.frame(WGCNA_res$net$colors),labels2colors(WGCNA_res$net$colors)) %>% tibble::rownames_to_column("gene")
  colnames(ME_genes)[2:3] <- c("cluster","Module")
  pos <- which(ME_genes$gene %in% c("ACP2","ADK","CDK10","DHFR","MAPKAPK3","MGST3","PKN3","THNSL2") )
  ME_genes[pos,]
  ME_green <- ME_genes$gene[ME_genes$Module == "green"]
  ME_green <- na.omit(convert(ME_green,fromtype = "SYMBOL",totype = "ENTREZID"))
  geneList <- t(WGCNA_res$dataExpr)[,1]
  geneList_names <- convert(names(geneList),fromtype = "SYMBOL",totype = "ENTREZID")
  pos <- which(!is.na(geneList_names))
  geneList <- geneList[pos]
  names(geneList) <- geneList_names[pos]
  ME_green_GOres <- FastGO(ME_green,geneList,organism = c("human","mouse")[1],
                           default.universe = F,classlevel = 2:2,
                           OrgDb  = NULL,  #c("org.Mm.eg.db","org.Hs.eg.db")
                           keyType = NULL,pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,qvalueCutoff  = 0.05,
                           cnet.showCategory = 5,verbose = TRUE,save.path = "ME_green_GO",names = "XieJing")
  
  ME_black <- ME_genes$gene[ME_genes$Module == "black"]
  ME_black <- na.omit(convert(ME_black,fromtype = "SYMBOL",totype = "ENTREZID"))
  ME_black_GOres <- FastGO(ME_black,geneList,organism = c("human","mouse")[1],
                           default.universe = F,classlevel = 2:2,
                           OrgDb  = NULL,  #c("org.Mm.eg.db","org.Hs.eg.db")
                           keyType = NULL,pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,qvalueCutoff  = 0.05,
                           cnet.showCategory = 5,verbose = TRUE,save.path = "ME_black_GO",names = "XieJing")
  
  ME_brown <- ME_genes$gene[ME_genes$Module == "brown"]
  ME_brown <- na.omit(convert(ME_brown,fromtype = "SYMBOL",totype = "ENTREZID"))
  ME_brown_GOres <- FastGO(ME_brown,geneList,organism = c("human","mouse")[1],
                           default.universe = F,classlevel = 2:2,
                           OrgDb  = NULL,  #c("org.Mm.eg.db","org.Hs.eg.db")
                           keyType = NULL,pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,qvalueCutoff  = 0.05,
                           cnet.showCategory = 5,verbose = TRUE,save.path = "ME_brown_GO",names = "XieJing")
  
  ME_turquoise <- ME_genes$gene[ME_genes$Module == "turquoise"]
  ME_turquoise <- na.omit(convert(ME_turquoise,fromtype = "SYMBOL",totype = "ENTREZID"))
  ME_turquoise_GOres <- FastGO(ME_turquoise,geneList,organism = c("human","mouse")[1],
                               default.universe = F,classlevel = 2:2,
                               OrgDb  = NULL,  #c("org.Mm.eg.db","org.Hs.eg.db")
                               keyType = NULL,pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,qvalueCutoff  = 0.05,
                               cnet.showCategory = 5,verbose = TRUE,save.path = "ME_turquoise_GO",names = "XieJing")
  save(ME_green_GOres,ME_black_GOres,ME_brown_GOres,ME_turquoise_GOres,file="D:/GSE5281/GO_res.rda")
}


