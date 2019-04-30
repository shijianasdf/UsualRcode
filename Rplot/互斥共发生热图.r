ME_h1<-function(M){
   pos<-which(colSums(M)!=0)
	 #M<-M[,pos]
	 flag<-paste(rep(letters,rep(10,length(letters))),0:9,sep="")
	 temp_c<-c()
    for(i in 1:dim(M)[2]){
         pos<-which(M[,i]==1)
         temp_c<-c(temp_c,paste(flag[pos],collapse=""))
    }
    names(temp_c)<-colnames(M)
    temp_c<-sort(temp_c)
    temp_sample<-names(temp_c)
    M<-M[,temp_sample]
	  return(M)
}

mat_hm_1<-ME_h1(cnv_matrix)
pdf(file="heatmap.pdf")
M<-as.matrix(mat_hm_1)
rownames(M)<-c("ANRIL","PTEN","PDGFRA","CHIC2")
colnames(M)<-NULL
library(gplots)
heatmap.2(M,trace="none",Colv=F,Rowv=F,
          col=c("white","red"),density.info="none",
          colsep=c(0,seq(0:ncol(M))), 
          rowsep=c(0,seq(0:nrow(M))), 
          sepcolor="gray", 
          sepwidth=c(0.0001,0.0001))
dev.off()

###
lnc<-"ENSG00000240498"
associated_gene<-rev(driver_lnc_all_related_genes$"ENSG00000240498")
common_sample<-intersect(colnames(cnv_M_lnc),colnames(cnv_mut_M))
cnv_matrix<-rbind(cnv_M_lnc[lnc,common_sample],cnv_mut_M[associated_gene,common_sample])