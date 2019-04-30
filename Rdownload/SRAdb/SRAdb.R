library(SRAdb)
sra_dbname <- file.path(system.file('extdata', package='SRAdb'), 'SRAmetadb_demo.sqlite')
sra_con <- dbConnect(dbDriver("SQLite"), sra_dbname)
## Get column descriptions
a <- colDescriptions(sra_con=sra_con)[1:5,]
## Convert SRA experiment accessions to other types
b <- sraConvert( in_acc=c(" SRR000137", "SRR000138 "), out_type=c('sample'), sra_con=sra_con )

SRXtoSRR.function <- function(SRXs){
    #@param SRXs
	#@param 返回列表(子元素是数据框)	
	library(SRAdb)
	sra_dbname <- '2017-06-05.SRAmetadb.sqlite'
	sra_con<- dbConnect(SQLite(), sra_dbname)
	SRRs <- vector(mode="list",length=length(SRXs))
	names(SRRs) <- SRXs
	for(i in 1:length(SRXs)){
		tempSRRs <- dbGetQuery(sra_con, paste("select library_layout,run_accession,experiment_accession from sra where experiment_accession = '",SRXs[i],"'",sep=""))
		SRRs[[i]] <- tempSRRs
	}
	return(SRRs)
}


