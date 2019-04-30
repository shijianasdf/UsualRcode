library(SRAdb)
sra_dbname <- file.path("/IData/CancerOMICS/BrainTumour(SRP145073)",'SRAmetadb.sqlite')
sra_con <- dbConnect(dbDriver("SQLite"), sra_dbname)
## Get column descriptions
a <- colDescriptions(sra_con=sra_con)[1:5,]
## Convert SRA experiment accessions to other types
b <- sraConvert( in_acc=c("SRP145073"), out_type=c('run'), sra_con=sra_con)
for(i in 1:nrow(b)){
	cat(paste("/pub5/xiaoyun/Software/aspera/connect/bin/ascp -i /pub5/xiaoyun/Software/aspera/connect/etc/asperaweb_id_dsa.openssh -QT -l300m",file.path("anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR",substring(b[i,2],1,6),b[i,2],paste(b[i,2],".sra",sep="")),"/IData/CancerOMICS/BrainTumour\\(SRP145073\\)/",sep=" "),"\n",file="/IData/CancerOMICS/BrainTumour(SRP145073)/download.sh",append=T);
}

## aspera:
/pub5/xiaoyun/Software/aspera/connect/bin/ascp -i /pub5/xiaoyun/Software/aspera/connect/etc asperaweb_id_dsa.openssh -QT -l300m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR713/SRR7138443/SRR7138443.sra /IData/CancerOMICS/BrainTumour\(SRP145073\)/ 



