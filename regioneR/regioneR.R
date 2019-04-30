library(regioneR)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
test <- data.frame(chr=c("chr1","chr2","chr3","chr4","chr1"),start=c(1,2,3,4,2),end=c(30,40,50,60,20),strand=c("+","-","-","+","-"))
test1 <- with(test,GRanges(chr,IRanges(start,end),strand))
getGenomeAndMask(Hsapiens, mask=NA)
randomizeRegions(test1, genome="hg19", mask=NULL, allow.overlaps=TRUE, per.chromosome=FALSE)

#Full genome sequences for Homo sapiens (Human) as provided by UCSC (hg19, Feb. 2009) and
#stored in Biostrings objects.
library(BSgenome.Hsapiens.UCSC.hg19)
seqnames(Hsapiens)
seqlengths(Hsapiens)
chr1 <- Hsapiens$chr1