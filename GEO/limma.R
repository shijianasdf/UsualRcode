#######################limma芯片差异##############################
library(limma)
design <- model.matrix( ~0 + factor( group_list ) )
colnames( design ) = levels( factor( group_list ) )
rownames( design ) = colnames( exprSet )
design
#制作比较矩阵"gbm-normal"意味着gbm相对normal求差异
contrast.matrix <- makeContrasts("gbm-normal",levels = design ) 
contrast.matrix
## design和contrast.matrix都是为了下面的差异分析制作的
DEGlimma <- function(exprSet,design,contrast.matrix){
  #library(limma)
  fit <- lmFit( exprSet, design )
  fit2 <- contrasts.fit( fit, contrast.matrix ) 
  fit2 <- eBayes( fit2 )
  tempOutput <- topTable( fit2, coef = 1, n = Inf )
  DEG <- na.omit(tempOutput)
  head(DEG)
  return(DEG)
}