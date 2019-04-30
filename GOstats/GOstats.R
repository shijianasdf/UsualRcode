#GO BP富集
#Before testing the terms whose children have already been tested, 
#we remove all genes annotated at significant children from the parent's gene list. 
#This continues until all terms have been tested
library(GOstats)
library(org.Hs.eg.db)
symbols <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2", 
       "ORM1",  "IGFBP1", "PTHLH", "GPC3",  "IGFBP3","TOB1",    "MITF",   "NDRG1", 
       "NR1H4", "FGFR3",  "PVR",   "IL6",   "PTPRM", "ERBB2",   "NID2",   "LAMB1", 
       "COMP",  "PLS3",   "MCAM",  "SPP1",  "LAMC1", "COL4A2",  "COL4A1", "MYOC",  
       "ANXA4", "TFPI2",  "CST6",  "SLPI",  "TIMP2", "CPM",     "GGT1",   "NNMT",
       "MAL",   "EEF1A2", "HGD",   "TCN2",  "CDA",   "PCCA",    "CRYM",   "PDXK",  
       "STC1",  "WARS",  "HMOX1", "FXYD2", "RBP4",   "SLC6A12", "KDELR3", "ITM2B");
selectedEntrezIds <- select(org.Hs.eg.db,keys=symbols,keytype="SYMBOL",columns="ENTREZID")
params <- new("GOHyperGParams",
               geneIds=selectedEntrezIds, #用户提交的基因集合
               universeGeneIds=entrezUniverse, #背景基因集合
               annotation="org.Hs.eg.db", #注释数据包
               ontology="BP",  #BP CC MF 富集到的GO term 一次只能选择一个
               pvalueCutoff=0.01,
               conditional=FALSE,
               testDirection="over") #over OR under
over <- hyperGTest(params) # return GOHyperGResult
GO <- summary(over)

#KEGG富集
params <- new("KEGGHyperGParams", 
                     geneIds=entrezIDs, 
                     universeGeneIds=universe, 
                     annotation="org.Hs.eg.db", 
                     categoryName="KEGG", 
                     pvalueCutoff=0.01,
                     testDirection="over")
over <- hyperGTest(params)
kegg <- summary(over)
