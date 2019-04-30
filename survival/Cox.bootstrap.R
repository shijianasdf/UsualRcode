#-------------------------------------
#'bootstrap cox回归模型的自证实
#-------------------------------------
mvcoxres_1 <- function(mod){
  hr <- round(summary(mod)$conf.int[, 1], 2)
  lower <- round(summary(mod)$conf.int[, 3], 2)
  upper <- round(summary(mod)$conf.int[, 4], 2)
  p <- round(summary(mod)$coefficients[, 5], 3)
  res <- data.frame(hr,lower,upper, p, stringsAsFactors = F)
  res$p[res$p < 0.001] <- "<.001"
  colnames(res) <- c("HR","lower","upper","p-value")
  return(res)
}
install.packages("bootStepAIC")
library(bootStepAIC)
library(survival)
library(help="bootStepAIC")
?boot.stepAIC
## 所有临床字段的结果
coxFit <- coxph(as.formula(paste("Surv(time,event)~",paste(colnames(COX.OS.all.gene)[-c(1:3)],collapse = "+"))), data=na.omit(COX.OS.all.gene))
aa <- boot.stepAIC(coxFit, data = na.omit(COX.OS.all.gene), B = 1000,direction = "both")
aaa <- boot.stepAIC(coxFit, data = na.omit(COX.OS.all.gene), B = 1000,direction = "backward")
a7 <- boot.stepAIC(coxFit, data = na.omit(COX.OS.all.gene), B = 1000,direction = "both")
a8 <- boot.stepAIC(coxFit, data = na.omit(COX.OS.all.gene), B = 1000,direction = "backward")
a10 <- boot.stepAIC(coxFit, data = na.omit(COX.OS.all.gene), B = 1000,direction = "both")
a11 <- boot.stepAIC(coxFit, data = na.omit(COX.OS.all.gene), B = 1000,direction = "both")
a10$Covariates 
a7$Covariates
aa$Covariates
a11$Covariates 
## 不加Residual临床字段的结果
coxFit1 <- coxph(as.formula(paste("Surv(time,event)~",paste(colnames(COX.OS.all.gene)[-c(1:3,10)],collapse = "+"))),data=na.omit(COX.OS.all.gene))
aaaa <- boot.stepAIC(coxFit1, data = na.omit(COX.OS.all.gene), B = 1000,direction = "backward")
aaaaa <- boot.stepAIC(coxFit1, data = na.omit(COX.OS.all.gene), B = 1000,direction = "both") 
a6 <- boot.stepAIC(coxFit1, data = na.omit(COX.OS.all.gene), B = 1000,direction = "both")
a9 <- boot.stepAIC(coxFit1, data = na.omit(COX.OS.all.gene), B = 1000,direction = "backward")
save(aa,aaa,a7,a8,aaaa,aaaaa,a6,a9,a10,a11,file="D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/bootstrap.RData")

load("D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/bootstrap.RData")
library(bootStepAIC)
bootstrap.resultTable.a7 <- lapply(a7$BootStepAIC,mvcoxres_1)
bootstrap.resultTable.a11 <- lapply(a11$BootStepAIC,mvcoxres_1)
bootstrap.resultTable.aa <- lapply(aa$BootStepAIC,mvcoxres_1)
bootstrap.resultTable.a10 <- lapply(a10$BootStepAIC,mvcoxres_1)
class(bootstrap.resultTable[[1]]) # data.frame
sj <- function(cox.resuts,rownam){
  ##cox.results
  ##rownam 
  result <- matrix(0,nrow=length(rownam),ncol=3,dimnames = list(rownam,c("mean HR","mean lower","mean upper")))
  for(j in 1:length(rownam)){
    temp <- list()
    for(i in 1:length(cox.resuts)){
      pos <- rownames(cox.resuts[[i]]) %in% rownam[j]
      if(any(pos)){
        temp[[i]] <- cox.resuts[[i]][pos,]
      }
    }
    tempMatrix <- do.call(rbind,temp)
    tempMatrix <- as.matrix(tempMatrix[,1:3])
    tempMatrix[tempMatrix>=7.261615e+24] <- Inf ## 将太大的值设为Inf
    tempMatrix[!is.finite(tempMatrix)] <- NA  ## 将Inf和-Inf改为NA
    tempRow <- colMeans(na.omit(tempMatrix))
    result[j,] <- tempRow
  }
  return(result)
}
rownam <- rownames(mvcoxres_1(coxFit))
boot.Mean.a7 <- sj(bootstrap.resultTable.a7,rownam) #a7
boot.Mean.a11 <- sj(bootstrap.resultTable.a11,rownam) #a11
boot.Mean.aa <- sj(bootstrap.resultTable.aa,rownam) #aa
boot.Mean.a10 <- sj(bootstrap.resultTable.a10,rownam) #a10
save(boot.Mean.a7,boot.Mean.a11,boot.Mean.aa,boot.Mean.a10,file="D:/R/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/bootstrap.mean.HR.RData")
  


## 利用ggforest画HR森林图
library(survminer)
ggforest(coxFit, data = na.omit(COX.OS.all.gene), 
         main = "Hazard ratio", 
         cpositions = c(0.10, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)
P <- ggforest(coxFit, data = na.omit(COX.OS.all.gene))
ggsave(P,filename = "D:/shijian.png")
