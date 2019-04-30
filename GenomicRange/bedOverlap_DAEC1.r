bed_overlap <- function(filepath,gene,disease){
	#用户上传文件路径
	#参数靶基因
	#参数变异
	#参数疾病
	#主表地址
	
	#tablepath <-"D:/TomcatZip/apache-tomcat-7.0.70/webapps/DiseaseEnhancer/RFunctions/diseaseEnh5.txt"
	tablepath <-"/pub5/xiaoyun/Software/tomcat7/webapps/DiseaseEnhancer/RFunctions/diseaseEnh5.txt"
	#读入数据(两个表格)
	#主表 chr start end TargetGene VarSite DiseaseType ...   TargetGene可能是,分隔的字符串   DiseaseType可能是;分隔的疾病
	mainTable <- read.csv(tablepath,sep="\t",header=T,stringsAsFactors=F)
	#flag <- (is.na(mainTable$chr) | is.na(mainTable$start) | is.na(mainTable$end))
	#mainTable <- mainTable[!flag,]
	#附表 chr start end strand
	queryTable  <- read.csv(filepath,sep="\t",header=F,stringsAsFactors=F)
	colnames(queryTable) <- c("chr", "start", "end")
	library(GenomicRanges)
	queryInput <- with(queryTable, GRanges(chr, IRanges(start,end))) #默认链信息为*
	#mainTable.Temp <- with(mainTable,GRanges(chr,IRanges(start,end)))
	mainTable.Temp <- with(mainTable,GRanges(chr,IRanges(start,end))) #默认链信息为*
	overlaps <- findOverlaps(queryInput, mainTable.Temp)
	#match_hit <- mainTable[names(mainTable.Temp)[subjectHits(overlaps)], c(1, 21, 2, 29:31, 10, 25)]
	match_hit <- mainTable[subjectHits(overlaps),]
	match_hit <- unique(match_hit)
	##############
	diseases <- unlist(strsplit(disease,", "))
	##############
	if(gene!="NULL"&disease!="NULL"){
	    match_hit <- match_hit[match_hit$diseasetype %in% diseases,] 
		match_hit <- match_hit[grep(gene,match_hit$targetgene),] 
	}
	if(gene!="NULL"&disease=="NULL"){
		match_hit <- match_hit[grep(gene,match_hit$targetgene),]
	}
	if(gene=="NULL"&disease!="NULL"){
	    match_hit <- match_hit[match_hit$diseasetype %in% diseases,]
	}
	if(gene=="NULL"&disease=="NULL"){
		match_hit <- match_hit
	}
	#if(gene!="NULL"&variant!="NULL"&disease!="NULL"){
	#    match_hit <- match_hit[match_hit$DiseaseType %in% diseases,] 
	#	match_hit <- match_hit[which(match_hit$TargetGene==gene & match_hit$varSite==variant ),] 
	#}
	#if(gene!="NULL"&variant=="NULL"&disease!="NULL"){
	#	match_hit <- match_hit[match_hit$DiseaseType %in% diseases,]
	#	match_hit <- match_hit[which(match_hit$TargetGene==gene),]
	#}
	#if(gene!="NULL"&variant!="NULL"&disease=="NULL"){
	#	match_hit <- match_hit[which(match_hit$varSite==variant & match_hit$TargetGene==gene),]
	#}
	#if(gene!="NULL"&variant=="NULL"&disease=="NULL"){
	#	match_hit <- match_hit[which(match_hit$TargetGene==gene),]
	#}
	#if(gene=="NULL"&variant!="NULL"&disease!="NULL"){
	#    match_hit <- match_hit[match_hit$DiseaseType %in% diseases,]
	#	match_hit <- match_hit[which(match_hit$varSite==variant),]
	#}
	#if(gene=="NULL"&variant=="NULL"&disease!="NULL"){
	#	match_hit <- match_hit[match_hit$DiseaseType %in% diseases,]
	#}
	#if(gene=="NULL"&variant!="NULL"&disease=="NULL"){
	#	match_hit <- match_hit[which(match_hit$varSite==variant),]
	#}
	#if(gene=="NULL"&variant=="NULL"&disease=="NULL"){
	#	match_hit <- match_hit
	#}
	#与java接口 需要as.matrix并且转置
	return(as.matrix(t(match_hit)))
}
