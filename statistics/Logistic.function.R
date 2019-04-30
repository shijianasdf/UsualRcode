#----------------------------------
  #'@author jianshi
  #'
  #'@description 利用给定的临床信息中的临床协变量（如：性别sex，年龄age等）
  #'进行单多因素以及逐步的logistic回归分析，输出分析结果。
  #'
  #'@param dat 即输入的临床信息
  #'注意：协变量若为离散型，必须转化为factor，且levels的第一项是参考项。
  #'@param y 即输入的响应变量
  #'@param step 即是否进行逐步logistic回归,目前step不好使，因为我不打算用
  #'@return 返回一个data.frame，包含单多因素logistic wald test的p值，风险比（OR），以及OR的
  #'置信区间（95% CI for OR）
  #'调用形式如下 My_logistic(test,"AJCC_NODES_PATHOLOGIC_PN",step=F) 
#----------------------------------
mvlogitres <- function(mod){
  or.ci <- paste0(
    round(exp(stats::coef(mod)[-1]), 2), " (",
    round(exp(stats::confint.default(mod)[-1, 1]), 2), " - ",
    round(exp(stats::confint.default(mod)[-1, 2]), 2), ")"
  )
  p <- round(summary(mod)$coefficients[-1, 4], 3)
  res <- data.frame(or.ci, p, stringsAsFactors = F)
  res$p[res$p < 0.001] <- "<.001"
  colnames(res) <- c("**OR (95% CI)**", "**p-value**")
  return(res)
}
My_logistic <- function(dat,y,step=c(T,F)){
  single.logistic <- function(dat,y){
    #STEP1:提取协变量名称
    covariates <- colnames(dat);
    covariates <- covariates[-match(y,covariates)]
    #STEP2:构建单因素分析的对象
    #glm(AJCC_NODES_PATHOLOGIC_PN~MSS_state,data=test)
    univ_formulas <- sapply(covariates,
                            function(x) paste('glm(',y,'~',x,',family = "binomial",data=dat)',sep=" "))
    #print(univ_formulas)
    #STEP3:单因素logistic分析
    univ_models <- lapply(univ_formulas, function(x){eval(parse(text=x))})

    #STEP4:提取有用信息
    univ_results <- lapply(univ_models,
                           function(x){
                             #提取p值，保留两位有效数字
                             #x <- glm(AJCC_METASTASIS_PATHOLOGIC_PM~APC,family = "binomial",data=na.omit(logist.clinical.M.data[[1]]))
                             #tmp <- summary(x)
                             print(tmp)
                             p.value <- round(tmp$coef[,4],4)
                             p.value[which(p.value < 0.0001)] <- "<0.0001"
                             tt <- exp(cbind(OR = coef(x), confint(x)))
                             variate <- rownames(tmp$coef)
                             #将所有值合并在一个矩阵中
                             all.data <- as.data.frame(cbind(variate, tt, p.value))[-1,]
                           })
    #STEP5: list转化为数据框输出
    univ.result <- do.call(rbind.data.frame, univ_results)
    colnames(univ.result) <- c('variate', 'univ OR','univ 2.5%','univ 97.5%','univ p value')
    rownames(univ.result) <- NULL
    return(univ.result)
  }
  multivariat.logistic <- function(dat,y,step){
    #STEP1:提取协变量名称
    covariates <- colnames(dat)
    covariates <- covariates[-match(y,covariates)]
    if(step==T){
      #STEP2:直接对所有选中的协变量进行逐步logistic
      command <- paste("step(glm(",y ,"~", paste(covariates, collapse="+"),",family = 'binomial',data=na.omit(dat)))", sep="")
      #print(command)
      mylogit <- eval(parse(text=command))
      mvlogitres(mylogit)
      tmp <- summary(mylogit) #得到最优logistic回归模型
      tt <- exp(cbind(OR = coef(mylogit), confint(mylogit)))
      coeff <- cbind.data.frame(rownames(tmp$coefficients),tmp$coefficients)
      colnames(coeff)[1] <- "variate" 
      variate <- rownames(tt)
      dd <- cbind.data.frame(variate, tt)
      mm <- merge(dd,coeff,by="variate",all.x=T)
      mm <- mm[,c(1,2,3,4,8)]
      mm[,5] <- round(mm[,5],digits = 4)
      mm[,5][which(mm[,5] < 0.0001)] <- "<0.0001"
      colnames(mm) <- c('variate', 'multi OR','multi 2.5%','multi 97.5%','multi p value')
      rownames(mm) <- NULL
      mm <- mm[-1,]
    }else{
      #STEP2:直接对所有选中的协变量进行多变量logistic
      command <- paste("glm(",y ,"~", paste(covariates, collapse="+"),",family = 'binomial', data=dat)", sep="")
      print(command)
      mylogit <- eval(parse(text=command))
      tmp <- summary(mylogit)
      print(tmp)
      tt <- exp(cbind(OR = coef(mylogit), confint(mylogit)))
      coeff <- cbind.data.frame(rownames(tmp$coefficients),tmp$coefficients)
      colnames(coeff)[1] <- "variate" 
      variate <- rownames(tt)
      dd <- cbind.data.frame(variate, tt)
      mm <- merge(dd,coeff,by="variate",all.x=T)
      mm <- mm[,c(1,2,3,4,8)]
      mm[,5] <- round(mm[,5],digits = 4)
      mm[,5][which(mm[,5] < 0.0001)] <- "<0.0001"
      colnames(mm) <- c('variate', 'multi OR','multi 2.5%','multi 97.5%','multi p value')
      rownames(mm) <- NULL
      mm <- mm[-1,]
    }
    return(mm)
  }
  uni.result <- single.logistic(dat,y)
  multi.result <- multivariat.logistic(dat,y,step)
  logistic.result <- merge(uni.result,multi.result,by="variate")
  return(logistic.result)
}
tt <- glm(AJCC_NODES_PATHOLOGIC_PN~SEX+MSS_state+AJCC_TUMOR_PATHOLOGIC_PT+AJCC_METASTASIS_PATHOLOGIC_PM+AJCC_PATHOLOGIC_TUMOR_STAGE+RESIDUAL_TUMOR+AGE+subtype+tumour_site_2+SMAD4,family = 'binomial', data=logist.clinical.N.data[[4]],control=list(maxit=1000))
mvlogitres(tt)
