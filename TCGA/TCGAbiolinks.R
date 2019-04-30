#https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

# In GDC database the clinical data can be retrieved from two sources:
#   
# indexed clinical: a refined clinical data that is created using the XML files.
# XML files
# There are two main differences:
#   
# XML has more information: radiation, drugs information, follow-ups, biospecimen, etc. So the indexed one is only a subset of the XML files
# The indexed data contains the updated data with the follow up informaiton. For example: if the patient is alive in the first time clinical data was collect and the in the next follow-up he is dead, the indexed data will show dead. The XML will have two fields, one for the first time saying he is alive (in the clinical part) and the follow-up saying he is dead. You can see this case here:
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)

library(help="TCGAbiolinks")
getGDCprojects()$project_id

#下面填入要下载的癌症种类
{
  request_cancer=c("PRAD","BLCA","KICH","KIRC","KIRP")
  for (i in request_cancer) {
    cancer_type=paste("TCGA",i,sep="-")
    print(cancer_type)
    #下载简单临床数据
    clinical <- GDCquery_clinic(project = cancer_type, type = "clinical")
    write.csv(clinical,file = paste(cancer_type,"clinical.csv",sep = "-"))
    
    #下载rna-seq的counts数据
    query <- GDCquery(project = cancer_type, 
                      data.category = "Transcriptome Profiling", 
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts")
    
    GDCdownload(query, method = "api", files.per.chunk = 100)
    expdat <- GDCprepare(query = query)
    count_matrix=assay(expdat)
    write.csv(count_matrix,file = paste(cancer_type,"Counts.csv",sep = "-"))
    
    #下载miRNA数据
    query <- GDCquery(project = cancer_type, 
                      data.category = "Transcriptome Profiling", 
                      data.type = "miRNA Expression Quantification", 
                      workflow.type = "BCGSC miRNA Profiling")
    
    GDCdownload(query, method = "api", files.per.chunk = 50)
    expdat <- GDCprepare(query = query)
    count_matrix=assay(expdat)
    write.csv(count_matrix,file = paste(cancer_type,"miRNA.csv",sep = "-"))
    
    #下载Copy Number Variation数据
    query <- GDCquery(project = cancer_type, 
                      data.category = "Copy Number Variation", 
                      data.type = "Copy Number Segment")
    
    GDCdownload(query, method = "api", files.per.chunk = 50)
    expdat <- GDCprepare(query = query)
    count_matrix=assay(expdat)
    write.csv(count_matrix,file = paste(cancer_type,"Copy-Number-Variation.csv",sep = "-"))
    
    #下载甲基化数据
    query.met <- GDCquery(project =cancer_type,
                          legacy = TRUE,
                          data.category = "DNA methylation")
    GDCdownload(query.met, method = "api", files.per.chunk = 300)
    expdat <- GDCprepare(query = query)
    count_matrix=assay(expdat)
    write.csv(count_matrix,file = paste(cancer_type,"methylation.csv",sep = "-"))
  }
}



#获取简单临床信息的方法 Get clinical indexed data
{
  cancer_type <- "TCGA-COAD"
  clinical <- GDCquery_clinic(project = cancer_type, type = "clinical")
  write.csv(clinical,file = paste(cancer_type,"clinical.csv",sep = "-"))
}
#Below are several examples fetching clinical data 
#directly from the clinical XML files.
{
  query <- GDCquery(project = "TCGA-READ", 
                    data.category = "Clinical", 
                    file.type = "xml"
                    #barcode = c("TCGA-RU-A8FL","TCGA-AA-3972"))
                    )
  GDCdownload(query,directory = "D:/GDCdata")
  #c("patient","admin","radiation","follow_up","drug","new_tumor_event")
  clinical <- GDCprepare_clinic(query, clinical.info = "patient",directory = "D:/GDCdata")
  write.csv(clinical,row.names = F,file="D:/TCGA_READ_patient.csv")
  clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug",directory = "D:/GDCdata")
  write.csv(clinical.drug,row.names = F,file="D:/TCGA_READ_clinical_drug.csv")
  clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation",directory = "D:/GDCdata")
  write.csv(clinical.radiation,row.names = F,file="D:/TCGA_READ_clinical_radiation.csv")
  clinical.admin <- GDCprepare_clinic(query, clinical.info = "admin",directory = "D:/GDCdata")
  write.csv(clinical.admin,row.names = F,file="D:/TCGA_READ_clinical_admin.csv")
  clinical.follow_up <- GDCprepare_clinic(query, clinical.info = "follow_up",directory = "D:/GDCdata")
  write.csv(clinical.follow_up,row.names = F,file="D:/TCGA_READ_clinical_follow_up.csv")
  clinical.stage_event <- GDCprepare_clinic(query, clinical.info = "stage_event",directory = "D:/GDCdata")
  write.csv(clinical.stage_event,row.names = F,file="D:/TCGA_READ_clinical_stage_event.csv")
  clinical.new_tumor_event <- GDCprepare_clinic(query, clinical.info = "new_tumor_event",directory = "D:/GDCdata")
  write.csv(clinical.new_tumor_event,row.names = F,file="D:/TCGA_READ_clinical_new_tumor_event.csv")
}
#下载TCGA-COAD的msi信息
#classifications: microsatellite-stable (MSS),
#low level MSI (MSI-L) if less than 40% of markers were altered 
#and high level MSI (MSI-H) if greater than 40% of markers were altered.
{
  query <- GDCquery(project = "TCGA-READ", 
                  data.category = "Other",
                  legacy = TRUE,
                  access = "open",
                  data.type = "Auxiliary test"
                 )  
  GDCdownload(query,directory = "D:/GDCdata")
  msi_results <- GDCprepare_clinic(query, "msi",directory = "D:/GDCdata")
  write.csv(msi_results,row.names = F,file="D:/TCGA_READ_msi_results.csv")
}
#Other useful code
#To get all the information for TGCA samples you can use the script below:
# This code will get all clinical indexed data from TCGA
{
  library(data.table)
  library(dplyr)
  library(regexPipes)
  clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>% 
    regexPipes::grep("TCGA",value=T) %>% 
    sort %>% 
    plyr::alply(1,GDCquery_clinic, .progress = "text") %>% 
    rbindlist
  readr::write_csv(clinical,path = paste0("all_clin_indexed.csv"))
  
  # This code will get all clinical XML data from TCGA
  getclinical <- function(proj){
    message(proj)
    while(1){
      result = tryCatch({
        query <- GDCquery(project = proj, data.category = "Clinical",file.type = "xml")
        GDCdownload(query)
        clinical <- GDCprepare_clinic(query, clinical.info = "patient")
        for(i in c("admin","radiation","follow_up","drug","new_tumor_event")){
          message(i)
          aux <- GDCprepare_clinic(query, clinical.info = i)
          if(is.null(aux) || nrow(aux) == 0) next
          # add suffix manually if it already exists
          replicated <- which(grep("bcr_patient_barcode",colnames(aux), value = T,invert = T) %in% colnames(clinical))
          colnames(aux)[replicated] <- paste0(colnames(aux)[replicated],".",i)
          if(!is.null(aux)) clinical <- merge(clinical,aux,by = "bcr_patient_barcode", all = TRUE)
        }
        readr::write_csv(clinical,path = paste0(proj,"_clinical_from_XML.csv")) # Save the clinical data into a csv file
        return(clinical)
      }, error = function(e) {
        message(paste0("Error clinical: ", proj))
      })
    }
  }
  clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>% 
    regexPipes::grep("TCGA",value=T) %>% sort %>% 
    plyr::alply(1,getclinical, .progress = "text") %>% 
    rbindlist(fill = TRUE) %>% setDF %>% subset(!duplicated(clinical))
  
  readr::write_csv(clinical,path = "all_clin_XML.csv")
}
# result: https://drive.google.com/open?id=0B0-8N2fjttG-WWxSVE5MSGpva1U
# Obs: this table has multiple lines for each patient, as the patient might have several followups, drug treatments,
# new tumor events etc...

# get subtype information
{
  dataSubt <- TCGAquery_subtype(tumor = "COAD")
  dataSubt1 <- TCGAquery_subtype(tumor = "READ")
  lgg.gbm.subtype <- TCGAquery_subtype(tumor = "lgg")
  subtypes <- PanCancerAtlas_subtypes()
  head(subtypes)
  table(subtypes$cancer.type)
  subtypes <- as.data.frame(subtypes)
  subtypes <- as.data.frame(subtypes[(subtypes$cancer.type == "COAD" | subtypes$cancer.type == "READ"),])
  write.csv(subtypes,file = "D:/coadREAD.subtype.csv",row.names = F)
}

#下载TCGA-COAD read-count
{
  query <- GDCquery(project = "TCGA-COAD", 
                    data.category = "Transcriptome Profiling", 
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts")
  
  GDCdownload(query, method = "api", files.per.chunk = 100,directory = "D:/R/Project/TCGA-COAD/ReadCount")
  expdat <- GDCprepare(query = query,directory = "D:/R/Project/TCGA-COAD/ReadCount")
  count_matrix=assay(expdat)
  write.csv(count_matrix,file = paste(cancer_type,"Counts.csv",sep = "-"))
}
?GDCquery
?GDCdownload
