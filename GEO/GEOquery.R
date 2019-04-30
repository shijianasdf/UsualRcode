## getGEO函数会下载GSE项目下的所有子集，得到的gset对象是一个list，‘GSE4290’只有一个项目，之后的实战会遇到多子集的情况
## ‘getGPL = T’会直接下载探针注释平台，如果报错，本文最后会附上，其他进行平台注释的方法
library("GEOquery")
GSE_name = 'GSE4290' #'GSE4290'
options('download.file.method.GEOquery' = 'libcurl') #windows系统
#下载表型，探针注释，以及处理好的表达矩阵
gset <- getGEO( GSE_name,getGPL = T,destdir="D:/R/Project/BioLearning/Transcription") 
#下载CEL原始数据
{
	getGEOSuppFiles(GSE_name, makeDirectory = TRUE, baseDir = "D:/R/Project/BioLearning/Transcription",fetch_files = TRUE, filter_regex = NULL)
	untar("D:/R/Project/BioLearning/Transcription/GSE31519/GSE31519_RAW.tar",exdir="D:/R/Project/BioLearning/Transcription/GSE31519")
}
save(gset, file = 'D:/R/Project/Learning/Transcriptiongset.Rdata')
options( stringsAsFactors = F )
gset = gset[[1]]
exprSet = exprs( gset ) ## “GEOquery”包中的exprs函数用来取出表达矩阵
pdata = pData( gset ) ## “GEOquery”包中的pData函数用来取出样本信息
group_list = as.character( pdata[, 35] )
dim( exprSet )
## [1] 54613   180
exprSet[ 1:3, 1:5 ]
##           GSM97793 GSM97794 GSM97795 GSM97796 GSM97797
## 1007_s_at  10178.1  10122.9   7826.6  11098.4   8668.9
## 1053_at      388.2    517.5    352.4    609.9    430.1

## 现在就应该对得到的矩阵有这样一个印象
## 这个矩阵有180列，列名是样本编号，54613行，行名是探针编号
## 矩阵的内容就是各个样本对应探针检测到的表达量
## 探针本身是没有意义的，所以我们下一步就要将探针转为基因名


## 筛选探针
GPL = gset@featureData@data ## 第三步getGEO函数下载数据时，直接下载了平台，GPL就是注释矩阵的平台数据
## 也就是探针和基因的对应关系

#加载下载的文件
gset <- getGEO(filename ="D:/R/Project/BioLearning/Transcription/GSE17537_series_matrix.txt.gz") 