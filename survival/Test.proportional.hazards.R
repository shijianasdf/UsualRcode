#'------------------------------------------------------------------------------------------------
#' @author jian shi
#' @description  Schoenfeld residuals to check the proportional hazards assumption
#' @param clinical.data 即输入的临床信息，至少需要包含四列信息Patient_ID表示临床样本ID；
#'        event表示生存事件结局；time表示生存事件的时间。默认对从第四列起的变量进行单因素分析，
#'        注意：协变量若为离散型，必须转化为factor，且levels的第一项是HR的参考项。
#' @example  
#' @return 图片和一个数据框
#'------------------------------------------------------------------------------------------------

Test.proportional.hazards <- function(clinical.data,picture.path){
  library("survival")
  library("survminer")
  #res.cox <- coxph(Surv(time, status) ~ age + sex + wt.loss, data =  lung)
  options(stringsAsFactors = FALSE);
  #STEP1:提取协变量名称
  covariates <- colnames(clinical.data)[-c(1:3)];
  #STEP2:直接对所有选中的协变量进行多因素分析
  command <- paste("coxph(Surv(time, event)~", paste(covariates, collapse="+"),", data=clinical.data)", sep="")
  print(command)
  res.cox <- eval(parse(text=command))
  test.ph <- cox.zph(res.cox)
  #pdf("D:/shijian.pdf")
  pdf(picture.path,height=2000,width = 2000)
  ggcoxzph(test.ph)
  dev.off()
  return(test.ph)
}

