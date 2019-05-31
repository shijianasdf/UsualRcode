#-----------------------------------------------------
#' @author :shijian
#' @依据driver基因的CCF值进行ROC曲线分析
#-----------------------------------------------------
## 导入数据
load("F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CRC.clinical.CCF.data.list.RData")

## 对临床数据表头进行清洗
head(CRC.clinical.CCF.data.list$APC)
Maxstate.byDriverGene.CCF <- lapply(CRC.clinical.CCF.data.list,function(x){
  x <- x[,c(1,10,11,18)] #对临床数据进行筛选
  colnames(x)[1:3] <- c("Patient_ID","event","time")
  #事件转换为数字类型
  x$event <- as.numeric(x$event)
  return(x) #如果不设置，默认返回最后一行
})
## 筛选KM生存显著的基因进行时间依赖ROC曲线分析
OS.significant <- c("ANK1","ARID1A","CASP8","GRIN2A","SMAD2")
Five.Maxstate.byDriverGene.CCF <- Maxstate.byDriverGene.CCF[names(Maxstate.byDriverGene.CCF) %in% OS.significant]
head(Five.Maxstate.byDriverGene.CCF[[1]])


TimeDeROC <- function(clinical.data,marker,cutoff){
  #clinical.data 临床数据
  #marker 
  #cutoff
  library(survivalROC)
  nobs <- nrow(clinical.data)
  ## 指标, METHOD = NNE
  Mayo4.1 <- survivalROC(Stime=clinical.data$time,  
                 status=clinical.data$event,      
                 marker = marker,     
                 predict.time = cutoff,
                 span = 0.25*nobs^(-0.20))
  pdf(paste("F:/Rsources/",names(clinical.data)[4],"_",cutoff,".pdf"))
  plot(Mayo4.1$FP, Mayo4.1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="red",   
       xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.1$AUC,3)), 
       ylab="TP",main=paste(names(clinical.data)[4],"at Month",cutoff,sep=" "),lwd=2,cex.main=1.3,cex.lab=1.2,cex.axis=1.2,font=1.2)
  abline(0,1)
  dev.off()
}
for(i in 1:5){
  TimeDeROC(Five.Maxstate.byDriverGene.CCF[[i]],Five.Maxstate.byDriverGene.CCF[[i]][,4],36)
}

## 时间依赖的ROC曲线分析R包一
library(survivalROC)
?survivalROC
nobs <- nrow(na.omit(Five.Maxstate.byDriverGene.CCF[[1]]))
cutoff <- 60
## MAYOSCORE 4, METHOD = NNE
temp <-  survivalROC(Stime = na.omit(Five.Maxstate.byDriverGene.CCF[[1]])$time, #生存时间 
                     status = na.omit(Five.Maxstate.byDriverGene.CCF[[1]])$event, #生存结局     
                     marker = na.omit(Five.Maxstate.byDriverGene.CCF[[1]])$SMAD2, #预测变量    
                     predict.time = cutoff, #时间点，选取3年(36个月)
                     span = 0.1,  
                     method = "KM")
plot(Mayo4.1$FP, Mayo4.1$TP, type="l", xlim=c(0,1), ylim=c(0,1),   
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo4.1$AUC,3)), 
     ylab="TP",main="Mayoscore 4, Method = NNE \n  Year = 3")
abline(0,1)

## 时间依赖的ROC曲线分析R包二
library(timeROC)
with(na.omit(Five.Maxstate.byDriverGene.CCF[[1]]),
     ROC <<- timeROC(T=time,#结局时间 
                     delta=event,#生存结局 
                     marker=SMAD2,#预测变量 
                     cause=1,#阳性结局赋值，比如死亡与否
                     weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
                     times=c(12,24,36,48,60),#时间点，选取12,24,36,48,60生存率
                     ROC = TRUE,
                     iid = TRUE)
)
## 画roc曲线
plot(ROC,time=36,col = "blue",add =FALSE)
## time是时间点，col是线条颜色，add指是否添加在上一张图中
plot(ROC,time=60,col = "red",add = F)
## AUC
ROC$AUC
## 置信区间
confint(ROC)

