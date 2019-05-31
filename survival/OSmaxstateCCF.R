#-----------------------------------------------------
#' @author :shijian
#' @依据driver基因的CCF值进行最大选择秩检验分析
#-----------------------------------------------------
## 读入每个基因打上CCF值标签的临床数据
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
## 筛选KM生存显著的基因进行最大选择秩检验
OS.significant <- c("ANK1","ARID1A","CASP8","GRIN2A","SMAD2")
Five.Maxstate.byDriverGene.CCF <- Maxstate.byDriverGene.CCF[names(Maxstate.byDriverGene.CCF) %in% OS.significant]

## 单因素maxstat
{
  library(maxstat)
  library(survival)
  mtHL.res <- list()
  for(i in names(Five.Maxstate.byDriverGene.CCF)){
    mtHL <- maxstat.test(eval(parse(text = paste("Surv(time, event) ~ ",i))),
                         data=na.omit(Five.Maxstate.byDriverGene.CCF[[i]]), smethod="LogRank", pmethod="Lau94",minprop=0,maxprop=0.999,abseps=0.01);
    mtHL.res[[i]] <- mtHL
  }
  save(mtHL.res,file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CRC.OS.ccf.maxstat.res.list.RData")
  load(file="F:/Rsources/Project/预后分析/克隆预后分析/结直肠癌克隆预后分析/Results/3.SurvivalKM/CRC.OS.ccf.maxstat.res.list.RData")
  table(na.omit(Five.Maxstate.byDriverGene.CCF[[2]])[,4])
  table(na.omit(Five.Maxstate.byDriverGene.CCF[[1]])[,4])
  table(na.omit(Five.Maxstate.byDriverGene.CCF[[3]])[,4])
  table(na.omit(Five.Maxstate.byDriverGene.CCF[[4]])[,4])
  table(na.omit(Five.Maxstate.byDriverGene.CCF[[5]])[,4])
}



## 将所有基因整合在一起做多因素maxstat所需要的数据
{
  library(dplyr)
  Five.Maxstate.byDriverGene.CCF.all <- Five.Maxstate.byDriverGene.CCF[[1]] %>% 
    full_join(Five.Maxstate.byDriverGene.CCF[[2]],by=c("Patient_ID","event","time")) %>% 
    full_join(Five.Maxstate.byDriverGene.CCF[[3]],by=c("Patient_ID","event","time")) %>% 
    full_join(Five.Maxstate.byDriverGene.CCF[[4]],by=c("Patient_ID","event","time")) %>% 
    full_join(Five.Maxstate.byDriverGene.CCF[[5]],by=c("Patient_ID","event","time"))
}
## 多因素maxstat,adjusted for more than one prognostic factor.
{
  library(maxstat)
  library(survival)
  multi.mstat <- maxstat.test(eval(parse(text=paste("Surv(time, event) ~",paste(colnames(Five.Maxstate.byDriverGene.CCF.all)[c(-1,-2,-3)],collapse = "+")))),
                              data=na.omit(Five.Maxstate.byDriverGene.CCF.all), smethod="LogRank", pmethod="Lau94",minprop=0,maxprop=0.999,abseps=0.01)
  plot(multi.mstat)
  
  opar <- par(no.readonly = T)
  par(mai=c(1.0196235, 1.0196235, 0.8196973, 0.4198450))
  multi.mstat <- maxstat.test(Surv(time, event) ~ SMAD2 + ARID1A + CASP8 + GRIN2A + ANK1,
                              data=na.omit(Five.Maxstate.byDriverGene.CCF.all), smethod="LogRank", pmethod="Lau94",abseps=0.01,minprop=0,maxprop=0.999)
  multi.mstat
  multi.mstat$p.value
  multi.mstat$maxstats
  multi.mstat$univp.values
  multi.mstat$maxstats[[1]]$estimate #cutoff
  multi.mstat$maxstats[[1]]$p.value #显著性P值
  plot(multi.mstat)
}
