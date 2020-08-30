## 基于TXDB外显子注释信息从bam文件获取基因的read count(利用summarizeOverlaps函数)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)
txdb <- makeTxDbFromGFF(file="D:/ThunderDownload/gencode.v35.annotation.gff3/gencode.v35.annotation.gff3", format="gff3", organism="Homo sapiens")
eByg <- exonsBy(txdb, by=c("gene")) # 获取TXDB基因外显子注释信息
bamFiles <- list.files("D:/Software/R-4.0.2/library/systemPipeRdata/extdata/bam",full.names = T,pattern="*.bam$") #获取bam文件地址
bfl <- BamFileList(bamFiles, yieldSize=50000, index=character())
counteByg <- summarizeOverlaps(eByg, bfl, mode="Union", 
                               ignore.strand=TRUE, 
                               # preprocess.reads=invertStrand,
                               inter.feature=FALSE, 
                               singleEnd=TRUE)
countDFeByg <- assays(counteByg)$counts
rpkmDFeByg <- apply(countDFeByg, 2, function(x) returnRPKM(counts=x, ranges=eByg))
write.table(countDFeByg, "results/countDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
write.table(rpkmDFeByg, "results/rpkmDFeByg.xls", col.names=NA, quote=FALSE, sep="\t")
