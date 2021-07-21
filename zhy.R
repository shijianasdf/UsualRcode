# 计算指定基因的在3类样本的差异表达情况
# 导入表达数据和突变数据以及表型数据
load("D:/Work/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.expression.data/colorectal.expression.ENSG.data.RData")
load("D:/Work/预后分析/克隆预后分析/结直肠癌克隆预后分析/colorectal.expression.data/colorectal.expression.data.RData")
load("D:/Work/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CRC.survival.data.list.RData")
COADREAD.ReadCount.matrix[1:4,1:4]
temp.COADREAD.ReadCount.matrix[1:4,1:4]
dim(temp.COADREAD.ReadCount.matrix)
dim(COADREAD.ReadCount.matrix)
CRC.survival.data.list$ANK1
genes <- c("IL4I1", "IDO1", "IFNG","MAPK12", "ANK1",
           "CASP8","SMAD2","ARID1A","CIITA","HLA-DMA","HLA-DPA1","HLA-DQA2","HLA-DOA","HLA-DPB1")

coadread.rsem <- read.table(file="D:/Work/预后分析/克隆预后分析/结直肠癌克隆预后分析/Data/colorectal.expression.data/COADREAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",sep="\t",
              header=T,fill=T,quote=NULL,stringsAsFactors=F)
coadread.rsem[1:4,1:4]
rownames(coadread.rsem) <- coadread.rsem[,1]
coadread.rsem <- coadread.rsem[,-1]
colnames(coadread.rsem) <- substring(colnames(coadread.rsem),1,12)
colnames(coadread.rsem) <- gsub("\\.","-",colnames(coadread.rsem))
coadread.rsem <- coadread.rsem[-1,]
library(stringr)
sub.coadread.rsem <- coadread.rsem[Fastmatch(genes,stringr::str_split(rownames(coadread.rsem),"\\|",simplify = T)[,1]),]
sub.coadread.rsem[1:4,1:4]
rownames(sub.coadread.rsem) <-  str_split(rownames(sub.coadread.rsem),"\\|",simplify = T)[,1]
common.patient <- intersect(colnames(sub.coadread.rsem),CRC.survival.data.list$ANK1$Patient_ID)
sub.coadread.rsem <- sub.coadread.rsem[,colnames(sub.coadread.rsem) %in% common.patient]
CRC.survival.data.list$ANK1 <- CRC.survival.data.list$ANK1[CRC.survival.data.list$ANK1$Patient_ID %in% common.patient,]
sub.coadread.rsem <- sub.coadread.rsem[,CRC.survival.data.list$ANK1$Patient_ID]
table(CRC.survival.data.list$ANK1$sample.label)

adata <- cbind.data.frame(t(sub.coadread.rsem),CRC.survival.data.list$ANK1$sample.label,colnames(sub.coadread.rsem))
adata$IL4I1 <- as.numeric(adata$IL4I1)
adata$IDO1 <- as.numeric(adata$IDO1)
adata$IFNG <- as.numeric(adata$IFNG)
adata$MAPK12 <- as.numeric(adata$MAPK12)
adata$ANK1 <- as.numeric(adata$ANK1)
adata$CASP8 <- as.numeric(adata$CASP8)
adata$SMAD2 <- as.numeric(adata$SMAD2)
adata$ARID1A <- as.numeric(adata$ARID1A)
adata$CIITA <- as.numeric(adata$CIITA)
adata$`HLA-DMA`<- as.numeric(adata$`HLA-DMA`)
adata$`HLA-DPA1` <- as.numeric(adata$`HLA-DPA1`)
adata$`HLA-DQA2` <- as.numeric(adata$`HLA-DQA2`)
adata$`HLA-DOA` <- as.numeric(adata$`HLA-DOA`)
adata$`HLA-DPB1` <- as.numeric(adata$`HLA-DPB1`)
colnames(adata)
colnames(adata)[15:16] <- c("clonality","sample")


library(ggpubr)
library(ggplot2)
library(ggsignif)
library(patchwork)
head(adata)
table(adata$clonality)
pos1 <- which(adata$clonality == "ANK1 WT")
pos2 <- which(adata$clonality == "ANK1 Clonal")
pos3 <- which(adata$clonality ==  "ANK1 Subclonal")

sts1 <- boxplot.stats(adata$IL4I1[pos1])$stats
sts2 <- boxplot.stats(adata$IL4I1[pos2])$stats
sts3 <- boxplot.stats(adata$IL4I1[pos3])$stats
sts <- range(c(sts1,sts2,sts3))
my_comparisons <- list( c("ANK1 WT", "ANK1 Subclonal"), c("ANK1 Subclonal", "ANK1 Clonal"), c("ANK1 WT", "ANK1 Clonal") )
p1 <- ggboxplot(adata, x = "clonality", y = "IL4I1",
          color = "clonality", palette = "jco",outlier.colour = NA) + 
          stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
          stat_compare_means(label.y = 50) #+    
         # coord_cartesian(ylim = c(sts[1]*0.95,sts[2]*1.05))    # Add global p-value 

p2 <- ggboxplot(adata, x = "clonality", y = "IDO1",
          color = "clonality", palette = "jco",outlier.colour = NA)+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

p3 <- ggboxplot(adata, x = "clonality", y = "IFNG",
          color = "clonality", palette = "jco",outlier.colour = NA)+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

p4 <- ggboxplot(adata, x = "clonality", y = "MAPK12",
          color = "clonality", palette = "jco",outlier.colour = NA)+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

p5 <- ggboxplot(adata, x = "clonality", y = "ANK1",
          color = "clonality", palette = "jco",outlier.colour = NA)+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

p6 <- ggboxplot(adata, x = "clonality", y = "CASP8",
                color = "clonality", palette = "jco",outlier.colour = NA)+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

p7 <- ggboxplot(adata, x = "clonality", y = "SMAD2",
                color = "clonality", palette = "jco",outlier.colour = NA)+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

p8 <- ggboxplot(adata, x = "clonality", y = "ARID1A",
                color = "clonality", palette = "jco",outlier.colour = NA)+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

p9 <- ggboxplot(adata, x = "clonality", y = "CIITA",
                color = "clonality", palette = "jco",outlier.colour = NA)+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

p10 <- ggboxplot(adata, x = "clonality", y = "HLA-DMA",
                color = "clonality", palette = "jco",outlier.colour = NA)+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

p11 <- ggboxplot(adata, x = "clonality", y = "HLA-DPA1",
                 color = "clonality", palette = "jco",outlier.colour = NA)+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

p12 <- ggboxplot(adata, x = "clonality", y = "HLA-DQA2",
                 color = "clonality", palette = "jco",outlier.colour = NA)+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

p13 <- ggboxplot(adata, x = "clonality", y = "HLA-DOA",
                 color = "clonality", palette = "jco",outlier.colour = NA)+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

p14 <- ggboxplot(adata, x = "clonality", y = "HLA-DPB1",
                 color = "clonality", palette = "jco",outlier.colour = NA)+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

p1 <- ggplot(data=adata,aes(x = clonality, y = IL4I1,color = clonality)) + 
  geom_boxplot(outlier.colour = NA) +
  theme_classic() +
  geom_signif(comparisons = my_comparisons, test = wilcox.test, step_increase = 0.2)

p2 <- ggplot(data=adata,aes(x = clonality, y = IDO1,color = clonality)) + 
  geom_boxplot(outlier.colour = NA) +
  theme_classic() +
  geom_signif(comparisons = my_comparisons, test = wilcox.test, step_increase = 0.2) #+

p3 <- ggplot(data=adata,aes(x = clonality, y = IFNG,color = clonality)) + 
  geom_boxplot(outlier.colour = NA) +
  theme_classic() +
  geom_signif(comparisons = my_comparisons, test = wilcox.test, step_increase = 0.2) #+



p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + p11 + p12 + p13 + p14
