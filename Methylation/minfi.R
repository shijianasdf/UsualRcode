if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("minfi")
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("minfiData")
library(minfi)
library(minfiData)
RGsetEx

## RGsetEx: RGChannelSet, 622,399 features
MsetEx <- preprocessRaw(RGsetEx)
## MsetEx: MethylSet, 485,512 features
GMsetEx <- mapToGenome(MsetEx)
## GMsetEx: GenomicMethylSet, 485,512 features

baseDir <- system.file("extdata", package = "minfiData")
list.files(baseDir)
list.files(file.path(baseDir, "5723646052"))

targets <- read.metharray.sheet(baseDir)
targets
sub(baseDir, "", targets$Basename)
RGset <- read.metharray.exp(targets = targets)
RGset
pd <- pData(RGset)
pd[,1:4]
RGset2 <- read.metharray.exp(file.path(baseDir, "5723646052"))
RGset3 <- read.metharray.exp(baseDir, recursive = TRUE)

targets2$Basename <- file.path(baseDir, targets2$Sentrix_ID, 
                               paste0(targets2$Sentrix_ID, 
                                      targets2$Sentrix_Position))

annotation(RGsetEx)

# 加载 minfi 包
library(minfi)
library(minfiData)
# 根据SampleSheet.csv 文件，读取所有样本的 .idat 文件
targets <- read.metharray.sheet(baseDir)
targets
sub(baseDir, "", targets$Basename)
rgSet <- read.metharray.exp(targets = targets)
# targets <- read.metharray.sheet("./", pattern="HumanMethylation450_Demo_Sample_Sheet.csv")
# rgSet <- read.metharray.exp(targets=targets)
# 计算探针的p值，过滤掉在任何以一个样本中p值大于0.01的探针
probeP <- detectionP(rgSet)
keep <-  apply(probeP, 1 , function(t){all(t < 0.01)})
rgSet <- rgSet[keep,]
# 过滤掉所有探针p值的均值大于0.05的样本
keep <- apply(probeP, 2, mean) < 0.05
rgSet <- rgSet[,keep]
# 预处理，背景降噪和归一化，注意，此处可以根据情况，替换成另外的算法
Gset.funnorm <- preprocessFunnorm(rgSet)
# 探针过滤，去除在性染色体上的探针
annotation <- getAnnotation(Gset.funnorm)
sex_probe <- rownames(annotation)[annotation$chr %in% c("chrX", "chrY")]
keep <-  !(featureNames(Gset.funnorm) %in% sex_probe)
Gset.funnorm <- Gset.funnorm[keep,]
# 去除覆盖了SNP位点的探针，如果感觉过滤掉的探针太多，可以适当上调SNP频率, 将maf的值变大 
GRset <- dropLociWithSnps(Gset.funnorm, snps=c("SBE","CpG"), maf=0)
# DMP, 探针水平的差异分析
beta  <- getBeta(GRset)
group <- pData(GRset)$Sample_Group
dmp   <- dmpFinder(beta, pheno = group  , type = "categorical")
head(dmp)
# DMR， 甲基化区域的差异分析
group <- pData(GRset)$Sample_Group
designMatrix <- model.matrix(~ group)
dmrs <- bumphunter(GRset,
                   design = designMatrix,
                   cutoff = 0.2, B=0, type="Beta")
head(dmrs$table[,1:4], n =3 )

