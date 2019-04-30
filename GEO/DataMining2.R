library(GEOquery)
load('D:/R/Project/BioLearning/Transcription/Glioma/gset.Rdata')
options( stringsAsFactors = F )
gset <- gset[[1]]
exprSet <- exprs( gset ) ## “GEOquery”包中的exprs函数用来取出表达矩阵
pdata <- pData( gset ) ## “GEOquery”包中的pData函数用来取出样本信息
group_list <- as.character( pdata[, 35] ) #得到样本的哪种病信息
table(group_list)
#astrocytoma, grade 2       astrocytoma, grade 3 
#4                          7                         19 
#glioblastoma, grade 4                  non-tumor oligodendroglioma, grade 2 
#77                         23                         38 
#oligodendroglioma, grade 3 
#12 
glioblastoma_exprSet <- exprSet[,grep("glioblastoma",group_list)] #恶性胶质瘤
nontumor_exprSet <- exprSet[,grep("non-tumor",group_list)] #正常样本
astrocytoma_exprSet <- exprSet[,grep("astrocytoma",group_list)]
oligodendroglioma_exprSet <- exprSet[,grep("oligodendroglioma",group_list)]

exprSet <- cbind(glioblastoma_exprSet,nontumor_exprSet)
group_list <- c(rep("case",length(grep("glioblastoma",group_list))),rep("control",length(grep("non-tumor",group_list)))) 

library(hgu133plus2.db)
id2symbol <- toTable(hgu133plus2SYMBOL)
exprSet <- exprSet[rownames(exprSet) %in% id2symbol[,1],]
id2symbol <- id2symbol[match(rownames(exprSet),id2symbol[,1]),]
temp <- by(exprSet,id2symbol[,2],function(x){return(x[which.max(rowMeans(x)),])}) #多个对应一个取max的probe
temp <- do.call(rbind,temp)
dim(temp)
exprSet <- temp
rownames(exprSet) <- id2symbol[,2][match(rownames(exprSet),id2symbol[,1])]
exprSet <- log2(exprSet+1)


library(limma)
design <- model.matrix(~0 + factor(group_list))#制作分组矩阵
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(exprSet)
design
contrast.matrix <- makeContrasts(paste(unique(group_list),sep="",collapse = "-"),levels = design )
contrast.matrix #制作比较矩阵"case-control"意味着case相对control
## design和contrast.matrix都是为了下面的差异分析制作的
DEGlimma <- function(exprSet,design,contrast.matrix){
  library(limma)
  fit <- lmFit( exprSet, design )
  fit2 <- contrasts.fit( fit, contrast.matrix ) 
  fit2 <- eBayes( fit2 )
  tempOutput <- topTable( fit2, coef = 1, n = Inf )
  DEG <- na.omit(tempOutput)
  head(DEG)
  return(DEG)
}
DEG <- DEGlimma(exprSet,design,contrast.matrix)
#为差异表达基因打标签
lab <- ifelse(abs(DEG$logFC) > 1 & DEG$adj.P.Val < 0.05,ifelse(DEG$logFC > 1,"up","down"),"no")
DEG <- cbind.data.frame(DEG,lab)

library(pheatmap)
choose_gene <- rownames(DEG)[DEG$lab %in% c("up","down")]
choose_matrix <- exprSet[choose_gene, ]
choose_matrix <- t(scale(t(choose_matrix)))
annotation_col <- data.frame(CellType = factor(group_list ))
rownames( annotation_col ) <- colnames( exprSet )
pheatmap( fontsize = 2, choose_matrix, annotation_col = annotation_col, show_rownames = F,annotation_legend = F);

library( "ggplot2" )
logFC_cutoff <- with(DEG, mean(abs(logFC ) ) + 2 * sd( abs( logFC ) ) )
logFC_cutoff
logFC_cutoff = 1 ## 文章中设置为1
this_tile <- paste0( 'Cutoff for logFC is ', round( logFC_cutoff, 3 ),                       '
                       The number of up gene is ', nrow(DEG[DEG$lab =='up', ] ),                      '
                       The number of down gene is ', nrow(DEG[DEG$lab =='down', ] ) )
volcano <- ggplot(data = DEG, mapping=aes( x = logFC, y = -log10(P.Value), color = lab)) +
    geom_point(alpha = 0.4, size = 1.75,mapping = aes(colour = lab)) +
    theme_set( theme_set( theme_bw( base_size = 15 ) ) ) +
    xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
    ggtitle( this_tile ) + theme( plot.title = element_text( size = 15, hjust = 0.5)) +
    scale_colour_manual(values = c('blue','black','red') )
print(volcano )


















 