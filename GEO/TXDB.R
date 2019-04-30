## load the library
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
## list the contents that are loaded into memory
ls('package:TxDb.Hsapiens.UCSC.hg38.knownGene')
help(package="TxDb.Hsapiens.UCSC.hg38.knownGene")
library(help="TxDb.Hsapiens.UCSC.hg38.knownGene")
class(TxDb.Hsapiens.UCSC.hg38.knownGene)
methods(class=class(TxDb.Hsapiens.UCSC.hg38.knownGene))
?TxDb
## show the db object that is loaded by calling it's name
TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txs <- transcriptsBy(TxDb, by="gene")
txs[["5594"]]
exby <- exonsBy(TxDb, by = "tx", use.names = TRUE)
getSeq(Hsapiens, txs[["5594"]])
getSeq(Hsapiens, exby[["uc002zvn.3"]])

TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
x <- isActiveSeq(TxDb)
x[seq_along(x)] <- FALSE
x["chrX"] <- TRUE
isActiveSeq(TxDb) <- x
utr5 <- fiveUTRsByTranscript(TxDb)
utr3 <- threeUTRsByTranscript(TxDb)
CDS <- cds(TxDb)
gene <- genes(TxDb)
exon <- exons(TxDb)
txs <- transcripts(TxDb)
utr5 <- unlist(utr5)
utr3 <- unlist(utr3)
prom <- promoters(gene, upstream=1000, downstream=0)
anno <- list("5UTR"=utr5,
             "3UTR"=utr3,
             "CDS" =CDS,
             "gene"=gene,
             "exon"=exon,
             "txs" =txs,
             "prom"=prom)
anno.rd <- lapply(anno, reduce)
cvglists.rd <- lapply(anno.rd, coverage)
d1 <- binOverRegions(cvglists.rd, TxDb)
d2 <- binOverGene(cvglists.rd, TxDb)
plotBinOverRegions(d1, main="binOverRegions.reduced")
plotBinOverRegions(d2, main="binOverGene.reduced")


