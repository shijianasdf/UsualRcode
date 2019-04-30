library(ArrayExpress)
EBI <- "E-MTAB-5885"
result <- queryAE(keywords = "cancer", species = "homo+sapiens")
mexp5885 <- getAE(EBI,path="D:/R/Project/BioLearning/Transcription/EBIfull/" ,type = "full") #下载全部数据
rawset <- ae2bioc(mageFiles = mexp5885)#converts local MAGE-TAB files into a AffyBatch 原始数据
pdata <- pData(rawset) #获取表型信息
#extracts the column names from processed MAGE-TAB and return them. The output is needed to call the function procset.
cn <- getcolproc(mexp5885) 
show(cn)
proset <- procset(mexp5885, cn[1]) #converts local MAGE-TAB files into an ExpressionSet.

#produces an AffyBatch, an ExpressionSet or a NChannelSet from a raw dataset from the ArrayExpress database.
#下载cel文件处理分析
ETABM25.affybatch <- ArrayExpress("E-MEXP-1416",path="D:/R/Project/BioLearning/Transcription/EBIraw/",save=T); #save设置T，下载信息会保存下来，否则加载进去后删除
print(ETABM25.affybatch)
#ExpressionFeatureSet (storageMode: lockedEnvironment)
#assayData: 1354896 features, 16 samples 
#  element names: exprs 
#protocolData
#  rowNames: 515G.CEL 459F.CEL ... 459A.CEL (16 total)
#  varLabels: exprs dates
#  varMetadata: labelDescription channel
#phenoData
#  rowNames: 515G.CEL 459F.CEL ... 459A.CEL (16 total)
#  varLabels: Source.Name Material.Type ... Factor.Value..sex. (35
#    total)
#  varMetadata: labelDescription channel
#featureData: none
#experimentData: use 'experimentData(object)'
#Annotation: pd.u133.x3p  
sampleNames(ETABM25.affybatch)
#[1] "515G.CEL" "459F.CEL" "515E.CEL" "515H.CEL" "515J.CEL" "515L.CEL"
# [7] "459D.CEL" "459C.CEL" "515K.CEL" "515D.CEL" "515F.CEL" "515B.CEL"
#[13] "515I.CEL" "515C.CEL" "459B.CEL" "459A.CEL"
colnames(pData(ETABM25.affybatch)) #表型信息sdrf
experimentData(ETABM25.affybatch) #实验信息idf
featureData(ETABM25.affybatch) #探针注释信息arf，但是提取不到
ETABM25.affybatch@annotation #[1] "pd.u133.x3p" 探针注释数据，在bioconductor下载
library(affy)
AEsetnorm <- rma(ETABM25.affybatch) #Background correcting Normalizing Calculating Expression
exprSet <- exprs(AEsetnorm) #提取表达