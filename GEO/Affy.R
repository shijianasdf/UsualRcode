library( "GEOquery" )
options( 'download.file.method.GEOquery' = 'libcurl' ) #windows系统
options( stringsAsFactors = F )
GSE_name = 'GSE31519' #disease
gset <- getGEO(GSE_name,getGPL = T,destdir="D:/R/Project/BioLearning/Transcription")
gset <- getGEO(filename="D:/R/Project/BioLearning/Transcription/GSE31519_series_matrix.txt.gz")
gset = gset[[1]]
exprSet = exprs( gset ) ## “GEOquery”包中的exprs函数用来取出表达矩阵
pdata = pData( gset ) ## “GEOquery”包中的pData函数用来取出样本信息
## 筛选探针
GPL = gset@featureData@data ## 第三步getGEO函数下载数据时，直接下载了平台，GPL就是注释矩阵的平台数据
## 也就是探针和基因的对应关系


#下载原始CEL文件
getGEOSuppFiles(GSE_name, makeDirectory = TRUE, baseDir = "D:/R/Project/BioLearning/Transcription",fetch_files = TRUE, filter_regex = NULL)
untar("D:/R/Project/BioLearning/Transcription/GSE31519/GSE31519_RAW.tar",exdir="D:/R/Project/BioLearning/Transcription/GSE31519")
library(affy)
affy.data <- ReadAffy(celfile.path="D:/R/Project/BioLearning/Transcription/GSE31519") #读入CEL文件
eset.rma <- rma(affy.data) #Background correcting Normalizing Calculating Expression
exprSet <- exprs(eset.rma) #提取表达

GSE_name1 = 'GSE20437'
gset1 <- getGEO(GSE_name1,getGPL = T,destdir="D:/R/Project/BioLearning/Transcription")
gset1 = gset1[[1]]
exprSet1 = exprs( gset1 ) ## “GEOquery”包中的exprs函数用来取出表达矩阵
pdata1 = pData( gset1 ) ## “GEOquery”包中的pData函数用来取出样本信息
## 筛选探针
GPL1 = gset1@featureData@data ## 第三步getGEO函数下载数据时，直接下载了平台，GPL就是注释矩阵的平台数据
## 也就是探针和基因的对应关系
#下载原始CEL文件
getGEOSuppFiles(GSE_name1, makeDirectory = TRUE, baseDir = "D:/R/Project/BioLearning/Transcription",fetch_files = TRUE, filter_regex = NULL)
untar("D:/R/Project/BioLearning/Transcription/GSE20437/GSE20437_RAW.tar",exdir="D:/R/Project/BioLearning/Transcription/GSE20437")
library(affy)
affy.data1 <- ReadAffy(celfile.path="D:/R/Project/BioLearning/Transcription/GSE20437")
eset.rma1 <- rma(affy.data1) #Background correcting Normalizing Calculating Expression
exprSet1 <- exprs(eset.rma1)


