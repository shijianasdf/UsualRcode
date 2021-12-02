library(foreach)
library(doParallel)
cores <- detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
re <- foreach(i=1:dim(rna_counts_matrix[pos,])[1], .combine="c") %dopar% {
  pos1 <- which(rna_counts_matrix[pos,-1][i,] > 0)
  pos2 <- which(atac_counts_matrix[overlap$Qindex,][i,] > 0)
  inter.pos <- intersect(pos1,pos2)
  ol <- length(inter.pos)
}
#stop cluster
stopCluster(cl)
