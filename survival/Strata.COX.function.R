#'-------------------------------------------------
#' @author jian shi
#' @description 分层COX分析
#' @param clinical.data 即输入的临床信息，至少需要包含四列信息Patient_ID表示临床样本ID；
#'        event表示生存事件结局；time表示生存事件的时间。默认对从第四列起的变量进行单因素分析，
#'        注意：协变量若为离散型，必须转化为factor，且levels的第一项是HR的参考项。
#' @param  clone.mutation.data 克隆突变数据
#' @example  Strata.Cox.function(CRC.survival.data.list[[1]],strata = "")
#'-------------------------------------------------
Strata.Cox.function <- function(clinical.data,strata){
  ###设置工作环境
  options(stringsAsFactors = FALSE);
  require(survival);
  
  #STEP1:提取协变量名称
  covariates <- colnames(clinical.data)[-c(1:3)];
  covariates <- covariates[-match(strata,covariates)];
  
  #STEP2:直接对所有选中的协变量进行分层多因素分析
  command <- paste("coxph(Surv(time, event)~", paste(covariates, collapse="+"),"+ strata(",strata,")",", data=clinical.data)", sep="");
  print(command);
  temp <- eval(parse(text=command)); 
  
  #STEP3:提取有用信息
  tmp <- summary(temp);
  
  #提取p值，保留两位有效数字
  p.value <- round(tmp$coefficients[ ,5], digits = 4);
  p.value[which(p.value < 0.0001)] <- "<0.0001"; 
  
  #提取beta值，这里得到的coefficients为矩阵，且有多行（每行对应一个协变量）
  #beta <- round(tmp$coefficients[ ,1], digits = 4);
  
  #提取风险比
  HR <- round(tmp$coefficients[ ,2], digits = 4);
  
  #提取95%置信区间上下界
  HR.confint.lower <- round(tmp$conf.int[ ,"lower .95"], digits = 4);
  HR.confint.upper <- round(tmp$conf.int[ ,"upper .95"], digits = 4);
  
  #合并风险比HR和置信区间为一个内容
  HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")");
  
  variate <- rownames(tmp$coefficients);
  
  #整合输出内容为data.frame
  multiv.result <- as.data.frame(cbind(variate, HR, p.value));
  colnames(multiv.result) <- c("variate", "strata multiv HR (95% CI for HR)", "strata multiv p value");
  rownames(multiv.result) <- NULL;
  
  return(multiv.result);	
}
clinical.data <- CRC.survival.data.list[[1]][,c(1,10,11,2,3,4,5,6,7,8,9,15,17,18)]
colnames(clinical.data)[c(2,3)] <- c("event","time")
clinical.data[,2] <- as.numeric(clinical.data[,2])
Strata.Cox.function(clinical.data,strata = "AJCC_PATHOLOGIC_TUMOR_STAGE")
coxph(Surv(time, event)~SEX+MSS_state+AJCC_TUMOR_PATHOLOGIC_PT+AJCC_NODES_PATHOLOGIC_PN+AJCC_METASTASIS_PATHOLOGIC_PM+RESIDUAL_TUMOR+AGE+subtype+tumour_site_2+sample.label+ strata(AJCC_PATHOLOGIC_TUMOR_STAGE), data=clinical.data)
