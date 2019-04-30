##getGEO函数会下载GSE项目下的所有子集，得到的gset对象是一个list，‘GSE4290’只有一个项目，之后的实战会遇到多子集的情况
## ‘getGPL = T’会直接下载探针注释平台，如果报错，本文最后会附上，其他进行平台注释的方法
library( "GEOquery" )
GSE_name <- 'GSE4290'
options( 'download.file.method.GEOquery' = 'libcurl' ) #windows系统
gset <- getGEO( GSE_name, getGPL = T ,destdir="D:/R/Project/BioLearning/Transcription")
save(gset, file = 'D:/R/Project/BioLearning/Transcription/Glioma/gset.Rdata')
options( stringsAsFactors = F )
gset <- gset[[1]]
exprSet <- exprs( gset ) ## “GEOquery”包中的exprs函数用来取出表达矩阵
pdata <- pData( gset ) ## “GEOquery”包中的pData函数用来取出样本信息
group_list <- as.character( pdata[, 35] ) #得到样本的哪种病信息
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
GPL <- gset@featureData@data ## 第三步getGEO函数下载数据时，直接下载了平台，GPL就是注释矩阵的平台数据
## 也就是探针和基因的对应关系

n_expr <- exprSet[ , grep( "non-tumor",         group_list )]
g_expr <- exprSet[ , grep( "glioblastoma",      group_list )]
a_expr <- exprSet[ , grep( "astrocytoma",       group_list )]
o_expr <- exprSet[ , grep( "oligodendroglioma", group_list )]
## 样本分组，新的表达矩阵只有normal和gbm样本
exprSet <- cbind( n_expr, g_expr )
group_list <- c(rep( 'normal', ncol( n_expr ) ), rep( 'gbm',    ncol( g_expr ) ) )
dim( exprSet )
exprSet[ 1:5, 1:5 ]
table( group_list )
save( exprSet, group_list, file = 'D:/R/Project/BioLearning/Transcription/Glioma/exprSet_by_group.Rdata')
## 筛选探针
GPL <- gset@featureData@data ## 第三步getGEO函数下载数据时，直接下载了平台，GPL就是注释矩阵的平台数据
## 也就是探针和基因的对应关系
colnames(GPL)
view( GPL )
## GPL的“ID”列是探针，‘Gene Symbol”列是基因名，将这两列取出来
ids <- GPL[ ,c( 1, 11 ) ]
## ‘Gene Symbol”列格式有些特殊
## 一个探针对应了多个基因
ids <- ids[which(!(ids$'Gene Symbol' == "")),]  #过滤掉没有基因注释的探针
a<-strsplit(as.character(ids[,2]), " /// ")
tmp <- mapply( cbind, ids[,1], a ) 
#ID2gene <- as.data.frame( tmp,stringsAsFactors =F)
ID2gene <- do.call(rbind.data.frame,tmp)
colnames( ID2gene ) <- c( "id", "gene" )
#第二步 得到geneID和gene类型的对应关系
#这个关系可以从人类的GTF文件中提取出来，GTF可以在这里下载：https://asia.ensembl.org/info/data/ftp/index.html
#下载后的文件经过下面的shell脚本处理，即可以得到基因与基因类型的对应关系
awk '{if(!NF || /^#/){next}}1' /public/reference/gtf/gencode/gencode.v25lift37.annotation.gtf | cut -f9 | sed 's/"//g'|sed 's/;//g' | awk '{ print $4"	"$8 }' |awk '{if(/^E/){next}}1'|  awk '{ print $2"	"$1 }' |sort -k 1 | uniq > gencode.v25lift37.annotation.gtf.gene2type
gene2type <- read.table( 'gencode.v25lift37.annotation.gtf.gene2type' )
colnames( gene2type ) = c( "gene", "type" )
dim( gene2type )
## [1] 58298     2
## 说明一共有五万多个基因的对应关系
head( gene2type )
##       gene           type
## 1  5S_rRNA           rRNA
## 2      7SK       misc_RNA
## 3 A1BG-AS1      antisense
sort( table( gene2type$type ) )
## lincRNA    processed_pseudogene    protein_coding
##    7601                   10197             20087
## 说明绝大多数的基因是编码蛋白的
## 剔除基因类型为“protein_coding”的对应关系
gene2type <- gene2type[ gene2type[,2] != 'protein_coding', ]
length( unique( gene2type$gene ) )
save( gene2type, file = 'Relationship_others_gene.Rdata' )

## 这一步根据ID2gene去除没有注释的探针
exprSet <- exprSet[ rownames(exprSet) %in% ID2gene[ , 1 ], ]
ID2gene <- ID2gene[ match(rownames(exprSet), ID2gene[ , 1 ] ), ]
dim( exprSet )
dim( ID2gene )
tail( sort( table( ID2gene[ , 2 ] ) ), n = 12L )
##相同基因的表达数据取最大值，五万多个探针，这一步相对会运行较长时间
MAX <- by( exprSet, ID2gene[ , 2 ],function(x) rownames(x)[ which.max( rowMeans(x) ) ] )
MAX <- as.character(MAX)
exprSet <- exprSet[ rownames(exprSet) %in% MAX , ]
rownames( exprSet ) <- ID2gene[ match( rownames( exprSet ), ID2gene[ , 1 ] ), 2 ]
exprSet <- log(exprSet) #log处理表达
save(exprSet, group_list, file = 'D:/R/Project/BioLearning/Transcription/Glioma/final_exprSet.Rdata')
##limma差异分析
library(limma)
design <- model.matrix( ~0 + factor( group_list ) )
colnames( design ) = levels( factor( group_list ) )
rownames( design ) = colnames( exprSet )
design       
#gbm normal
#GSM97812   0      1
#GSM97820   0      1
#GSM97825   0      1
#GSM97827   0      1
#GSM97833   0      1
#GSM97840   0      1
#GSM97848   0      1
contrast.matrix <- makeContrasts("gbm-normal",levels = design )
contrast.matrix
## design和contrast.matrix都是为了下面的差异分析制作的
fit <- lmFit( exprSet, design )
fit2 <- contrasts.fit( fit, contrast.matrix ) 
fit2 <- eBayes( fit2 )
nrDEG <- topTable( fit2, coef = 1, n = Inf ) 
#表达矩阵热图
library( "pheatmap" )
choose_gene <- head( rownames( nrDEG ), 50 )
choose_matrix <- exprSet[ choose_gene, ]
choose_matrix <- t( scale( t( choose_matrix ) ) )
annotation_col <- data.frame( CellType = factor( group_list ) )
rownames( annotation_col ) <- colnames( exprSet )
pheatmap( fontsize = 2, choose_matrix, annotation_col = annotation_col, show_rownames = F,      annotation_legend = F);
#差异火山图
library( "ggplot2" )
logFC_cutoff <- with( nrDEG, mean( abs( logFC ) ) + 2 * sd( abs( logFC ) ) )
logFC_cutoff
logFC_cutoff = 1 ## 文章中设置为1
{
  nrDEG$change = as.factor(ifelse( nrDEG$P.Value < 0.05 & abs(nrDEG$logFC) > logFC_cutoff,
                                  ifelse( nrDEG$logFC > logFC_cutoff , 'UP', 'DOWN' ),'NOT'))
  save( nrDEG, file = "D:/R/Project/BioLearning/Transcription/nrDEG.Rdata" )
  this_tile <- paste0( 'Cutoff for logFC is ', round( logFC_cutoff, 3 ),                       '
The number of up gene is ', nrow(nrDEG[ nrDEG$change =='UP', ] ),                      '
The number of down gene is ', nrow(nrDEG[ nrDEG$change =='DOWN', ] ) )
  volcano = ggplot(data = nrDEG, mapping=aes( x = logFC, y = -log10(P.Value), color = change)) +
                   geom_point(alpha = 0.4, size = 1.75) +
                   theme_set( theme_set( theme_bw( base_size = 15 ) ) ) +
                   xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
                   ggtitle( this_tile ) + theme( plot.title = element_text( size = 15, hjust = 0.5)) +
                   scale_colour_manual(values = c('blue','black','red') )
  print(volcano )
  ggsave(volcano, filename = 'D:/R/Project/BioLearning/Transcription/Glioma/volcano.png' )
  dev.off()
}





















































