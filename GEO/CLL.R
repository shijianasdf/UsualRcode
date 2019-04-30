library(CLL)
library(help="CLL")
data(sCLLex)
sCLLex
?sCLLex
exprSet <- exprs(sCLLex)
group_list <- pData(sCLLex)$Disease

library(hgu95av2.db)
library(help="hgu95av2.db")
ID2gene <- toTable(hgu95av2SYMBOL)
##相同基因的表达数据取最大值，五万多个探针，这一步相对会运行较长时间
exprSet <- exprSet[ rownames(exprSet) %in% ID2gene[ , 1 ], ]
ID2gene <- ID2gene[ match(rownames(exprSet), ID2gene[ , 1 ] ), ]
MAX <- by( exprSet, ID2gene[ , 2 ],function(x) rownames(x)[ which.max( rowMeans(x) ) ] )
MAX <- as.character(MAX)
exprSet <- exprSet[ rownames(exprSet) %in% MAX , ]
rownames( exprSet ) <- ID2gene[ match( rownames( exprSet ), ID2gene[ , 1 ] ), 2 ]

##limma差异分析
library(limma)
design <- model.matrix( ~0 + factor( group_list ) )
colnames( design ) <- levels( factor( group_list ) )
rownames( design ) <- colnames( exprSet )
design  
contrast.matrix <- makeContrasts("progres.-stable",levels = design )
contrast.matrix
## design和contrast.matrix都是为了下面的差异分析制作的
fit <- lmFit( exprSet, design )
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)
DEG <- topTable(fit2, coef = "progres.-stable", n = Inf)
DEG1 <- topTable(fit2, coef = 1, n = Inf) #两个方法是一样的

#表达矩阵热图
library( "pheatmap" )
choose_gene <- head( rownames(DEG), 100 )
choose_matrix <- exprSet[ choose_gene, ]
choose_matrix <- t( scale( t( choose_matrix ) ) )
annotation_col <- data.frame( CellType = factor( group_list ) )
rownames( annotation_col ) <- colnames( exprSet )
pheatmap(choose_matrix, annotation_col = annotation_col, show_rownames = F,annotation_legend = F,fontsize = 2);

#差异火山图
library( "ggplot2" )
logFC_cutoff <- with(DEG, mean(abs(logFC)) + 2 * sd(abs(logFC)))
logFC_cutoff
logFC_cutoff <- 1 ## 文章中设置为1
{
    #为差异表达基因打标签
    DEG$lab <- as.factor(ifelse(abs(DEG$logFC) > logFC_cutoff & DEG$adj.P.Val < 0.05,ifelse(DEG$logFC > logFC_cutoff,"UP","DOWN"),"NO"))
    this_tile <- paste0('Cutoff for logFC is ', round(logFC_cutoff, 3 ),'
    The number of up gene is ', nrow(DEG[ DEG$lab =='UP', ] ),'
    The number of down gene is ', nrow(DEG[ DEG$lab =='DOWN', ] ) )
    volcano <- ggplot(data = DEG, mapping=aes(x = logFC, y = -log10(adj.P.Val), color = lab)) +
    geom_point(alpha = 0.4, size = 1.75) +
    theme_set( theme_set( theme_bw( base_size = 15 ) ) ) +
    xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
    ggtitle( this_tile ) + theme( plot.title = element_text( size = 15, hjust = 0.5)) +
    scale_colour_manual(values = c('blue','black','red') )
    print(volcano)
}

##加载原始CEL文件进行分析
data(CLLbatch)
?CLLbatch
library(affy)
exprSet.cel <- rma(CLLbatch)
exprSet.cel <- exprs(exprSet.cel)
dim(exprSet.cel)
colnames(exprSet.cel)
#样本临床信息
data(disease)
disease
group_list.cel <- na.omit(disease$Disease)

library(hgu95av2.db)
ID2gene <- toTable(hgu95av2SYMBOL)
##相同基因的表达数据取最大值，五万多个探针，这一步相对会运行较长时间
exprSet.cel <- exprSet.cel[rownames(exprSet.cel) %in% ID2gene[ , 1 ], ]
ID2gene <- ID2gene[match(rownames(exprSet.cel), ID2gene[ , 1 ] ), ]
MAX <- by(exprSet.cel,ID2gene[ , 2 ],function(x)rownames(x)[which.max(rowMeans(x))])
MAX <- as.character(MAX)
exprSet.cel <- exprSet.cel[ rownames(exprSet.cel) %in% MAX , ]
rownames(exprSet.cel) <- ID2gene[match( rownames( exprSet.cel), ID2gene[, 1 ]), 2]
exprSet.cel <- exprSet.cel[,-1]
##limma差异分析
library(limma)
design <- model.matrix( ~0 + factor( group_list.cel ) )
colnames( design ) <- levels( factor( group_list.cel ) )
rownames( design ) <- colnames( exprSet.cel )
design  
contrast.matrix <- makeContrasts("progres.-stable",levels = design )
contrast.matrix
## design和contrast.matrix都是为了下面的差异分析制作的
fit <- lmFit( exprSet.cel, design )
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)
DEG <- topTable(fit2, coef = 1, n = Inf)

#差异火山图
library( "ggplot2" )
logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2 * sd(abs(logFC)))
logFC_cutoff
logFC_cutoff <- 1 ## 文章中设置为1
{
    #为差异表达基因打标签
    DEG$lab <- as.factor(ifelse(abs(DEG$logFC) > logFC_cutoff & DEG$adj.P.Val < 0.05,ifelse(DEG$logFC > logFC_cutoff,"UP","DOWN"),"NO"))
    this_tile <- paste0('Cutoff for logFC is ', round(logFC_cutoff, 3 ),'
    The number of up gene is ', nrow(DEG[ DEG$lab =='UP', ] ),'
    The number of down gene is ', nrow(DEG[ DEG$lab =='DOWN', ] ) )
    volcano <- ggplot(data = DEG, mapping=aes(x = logFC, y = -log10(adj.P.Val), color = lab)) +
    geom_point(alpha = 0.4, size = 1.75) +
    theme_set( theme_set( theme_bw( base_size = 15 ) ) ) +
    xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
    ggtitle( this_tile ) + theme( plot.title = element_text( size = 15, hjust = 0.5)) +
    scale_colour_manual(values = c('blue','black','red') )
    print(volcano)
}













