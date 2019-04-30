library(airway)
data(airway)
airway
##RNAseq readcount表达信息64102*8
exprSet <- assay(airway) 
##临床信息
sampleInfo <- colData(airway) 
##分组信息
group_list <- sampleInfo$dex #untrt trt   untrt trt   untrt trt   untrt trt
##过滤掉少于50%样本表达的基因
temp <- apply(exprSet,1,function(x){
  t <- x > 0
  ifelse(sum(t)>0.5*length(x),T,F) 
})
exprSet <- exprSet[temp,] #42735*8
######################################################################
###################      Firstly for DEseq2      #####################
######################################################################
library(DESeq2)
(colData <- data.frame(row.names=colnames(exprSet), group_list=group_list))
colData
#group_list
#SRR1039508      untrt
# SRR1039509        trt
# SRR1039512      untrt
# SRR1039513        trt
# SRR1039516      untrt
# SRR1039517        trt
# SRR1039520      untrt
# SRR1039521        trt
group_list
# untrt trt   untrt trt   untrt trt   untrt trt
#产生DESeq对象
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)
#计算差异
dds <- DESeq(dds)
#替换异常值
dds <- replaceOutliersWithTrimmedMean(dds)
#提取结果
#提取结果trt相对于untrt的差异表达
res <- results(dds, contrast=c("group_list","trt","untrt")) 
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG=as.data.frame(resOrdered)
#提取结果untrt相对于trt的差异表达
res1 <- results(dds, contrast=c("group_list","untrt","trt"))
resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)
DEG1=as.data.frame(resOrdered1)


######################################################################
###################      Then  for edgeR        #####################
######################################################################
library(edgeR)
dge <- DGEList(counts=exprSet,group=factor(group_list))
dge <- calcNormFactors(dge)

design <- model.matrix(~0+factor(group_list))
rownames(design)<-colnames(dge)
colnames(design)<-levels(factor(group_list))
design

dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit,  contrast=c(1,-1))
nrDEG <- topTags(lrt, n=nrow(exprSet))
nrDEG <- as.data.frame(nrDEG)
lrt1 <- glmLRT(fit,  contrast=c(-1,1))
nrDEG1 <- topTags(lrt1, n=nrow(exprSet))
nrDEG1 <- as.data.frame(nrDEG1)


######################################################################
###################      Then  for limma/voom        #################
######################################################################
suppressMessages(library(limma))
library(edgeR)
design <- model.matrix(~0+factor(group_list))
colnames(design)<-levels(factor(group_list))
rownames(design)<-colnames(exprSet)
dge <- DGEList(counts=exprSet)
dge <- calcNormFactors(dge)
v <- voom(dge,design)
fit <- lmFit(v, design)
group_list
cont.matrix=makeContrasts(contrasts='trt-untrt',levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)
tempOutput = topTable(fit2, coef='trt-untrt', n=Inf)
DEG_limma_voom = na.omit(tempOutput)
