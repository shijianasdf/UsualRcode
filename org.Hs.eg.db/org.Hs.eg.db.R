library(org.Hs.eg.db)
#列举所有ENSEMBLTRANS
keys(org.Hs.eg.db,keytype="ENSEMBLTRANS")
#列举所有keytypes
keytypes(org.Hs.eg.db)
#ENSEMBL转换为gene symbol和gene name
ensids <-c("ENSG00000144644","ENSG00000159307","ENSG00000144485")
cols <- c("SYMBOL","GENENAME")
select(org.Hs.eg.db,keys=ensids,keytype = "ENSEMBL",columns = cols)
#gene symbol转换为其它
select(org.Hs.eg.db,keys="BRCA1",keytype="SYMBOL",columns=c("ENSEMBL","UNIGENE","ENTREZID","CHR","GO","GENENAME"))
mget("BRCA1", org.Hs.egSYMBOL2EG, ifnotfound=NA);