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
intersect(WGCNA_res$Hub$Block1$HubGenes,c("ACP2","ADK","CDK10","DHFR","MAPKAPK3","MGST3","PKN3","THNSL2")) 
intersect(colnames(WGCNA_res$dataExpr),c("ACP2","ADK","CDK10","DHFR","MAPKAPK3","MGST3","PKN3","THNSL2")) 

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
