####################################################################################
# 整体功能简介：
# 计算预后模型比较指标：C-index，AUC，IDI和NRI，AIC的值；同时在需要的情况下返回指标相关
# 的95%CI和比较的p值
###################################################################################
# 重要使用函数简介
#（1）cIndex：计算C-index的值，95%CI和C-index比较的p值
#（2）survivalROCauc：计算AUC的值，95%CI和比较的p值，同时绘制ROC曲线
#（3）NRI.IDI：计算NRI和IDI的值，95%CI和p值
#（4）aic：计算AIC的值
###################################################################################
#  论坛网址：

###################################################################################
# 创建作者: 赵二杰
# 日期：2019/01/07
###################################################################################


###################################################################################
# 修改日志
###################################################################################
# 1) 日期：              修改者：
#    简要修改内容：
# 
# 
###################################################################################
# 2) 日期：               修改者：
#    简要修改内容：
# 
# 
###################################################################################


##---------------------------------------------------------------------------------------------
#  calculate the C-index for a cox model
#'@description 计算生存模型的C-index或者两个模型比较的p值
#'@param clinical.data: 临床变量信息，至少包含model中的变量信息
#'@param time: 数值型向量，患者对应的生存时间
#'@param event: 患者对应的生存状态，通常0=alive, 1=dead
#'@param models：列表，每一个元素是一个字符型向量，包含一个model所有变量对应的列名,至少包含两个元素
#'@param diff: 逻辑值，diff = FALSE不计算两个模型比较的p值，只返回对应的C-index值；diff = TRUE返回模型对应的C-index值及模型比较的p值
#'@return 返回一个data.frame，包含C-index值，置信区间，比较的p值（当diff = T时）
cIndex <- function(clinical.data, time, event, models, diff = TRUE){
  
  options(stringsAsFactors = FALSE)

  # load survival, Hmisc package
  suppressPackageStartupMessages(require(survival))
  suppressPackageStartupMessages(require(Hmisc))

  ##-------------------------------
  # CoxphObject：create a coxph object
  CoxphObject <- function(data, time, event, variables){
    if (length(variables) == 1) {
	  formula <- as.formula(paste("Surv(time, event)~", variables))
	} else {
	  formula <- as.formula(paste("Surv(time, event)~", paste(variables, collapse=" + ")))
	}      
    coxph.object <-coxph(formula, data = data)

    return(coxph.object)
  }

  ##-------------------------------
  # 输出的Cindex和置信区间
  c.ci <- function(rcorrobj){
    CIndex <- round(rcorrobj['C Index'], digits = 4)
    se     <- rcorrobj['S.D.']/2
    Lower <- round(CIndex - 1.96*se, digits = 4)
    Upper <- round(CIndex + 1.96*se, digits = 4)
    result <- c(CIndex, Lower, Upper)
    names(result) <- c("C-Index", "Lower", "Upper")

    return(result)
  }
  
  # 计算每个模型的C-index值和置信区间
  coxph.list <- lapply(models, function(x){CoxphObject(data = clinical.data, time = time, event = event, variables = x)})
  pred.models.coxph <- lapply(coxph.list, function(x){predict(x, type = "risk")})
  models.result <- lapply(pred.models.coxph, function(x){rcorr.cens(-x, Surv(time = time, event = event))})
  
  models.filter.result <- lapply(models.result, function(x){c.ci(rcorrobj = x)})
  
  # 是否进行C-index1的比较
  if (diff == FALSE) {
    result <- do.call(rbind, models.filter.result)
    conf.interval <- paste(result[, 'Lower'], result[, 'Upper'], sep = "-")
    result <- data.frame(result[, 'C-Index'], conf.interval)
    colnames(result) <- c('C-index', '95%CI')

    return(result)
  } else {
    # 计算比较的p值，都是和第一个模型比较
    compare.cindex <- lapply(pred.models.coxph[-1], function(x){rcorrp.cens(pred.models.coxph[[1]], x, Surv(time = time, event = event))})
    p.value <- c("-", unlist(lapply(compare.cindex, function(x)(round(2*(1 - pnorm(abs(x[['Dxy']] / x[['S.D.']]))), digits=4))))) 

    result <- do.call(rbind, models.filter.result)
    conf.interval <- paste(result[, 'Lower'], result[, 'Upper'], sep = "-")
    result <- data.frame(result[, 'C-Index'], conf.interval, p.value)
    colnames(result) <- c('C-index', '95%CI', 'p_value')

    return(result)
  }
}  



#-------------------------------------------------------------------
# 计算生存模型的AUC和95%CI，以及模型比较时对应的p值
#'@description 绘制cox模型对应的ROC曲线，计算对应的AUC值
#'@param clinical.data: 临床变量信息，至少包含model中的变量信息
#'@param time:数值型向量，患者对应的生存时间
#'@param event:数值型向量，患者对应的生存状态，注意必须是0=alive, 1=dead
#'@param models：列表，每一个元素是一个字符型向量，包含一个model所有变量对应的列名,至少包含两个元素
#'@param upper.time：数值，绘制ROC曲线时截止的时间点，如：1年，3年。注意：数值必须与time变量的单位一致
#'@param col: 字符型向量，表示每个模型对应曲线的颜色，当col = NULL(默认)时，随机生成颜色
#'@param pdf.file：对应图像输出路径，默认为NULL
#'@return 返回两个信息：1,ROC曲线图像；2,返回一个data.frame，包含每个模型对应的AUC和95%CI，以及auc比较的p值
survivalROCauc <- function(clinical.data, time, event, models, upper.time, col = NULL, pdf.file = NULL){

  # load survival, pROC package
  suppressPackageStartupMessages(require(survival))
  suppressPackageStartupMessages(require(pROC))

  ##-------------------------------
  # CoxphObject：create a coxph object
  CoxphObject <- function(data, time, event, variables){
    if (length(variables) == 1) {
	  formula <- as.formula(paste("Surv(time, event)~", variables))
	} else {
	  formula <- as.formula(paste("Surv(time, event)~", paste(variables, collapse=" + ")))
	}      
    coxph.object <-coxph(formula, data = data)

    return(coxph.object)
  }
  
  models.list <- lapply(models, function(x){CoxphObject(data = clinical.data, time = time, event = event, variables = x)})

  # Filter survival data  ?为什么过滤,活着的并且时间小于upper.time
  index <- which(time < upper.time & event == 0) 
  survival.data <- clinical.data[-index, ]

  # predict risk value of models
  risk.list <- lapply(models.list, function(x){predict(x, newdata = survival.data, type = 'risk')})
  
  #Label the sample, dead in upper.time is 1, alive in upper.time is 0
  label <- as.numeric(time[-index] <= upper.time) 

  rocobj.list <- lapply(risk.list, function(x){roc(label ~ x, ci = T, auc = T)})
  # Calculate the auc value for each model
  auc.list <- lapply(rocobj.list, function(x){auc(x)[1]})

  # Calculate the 95% CI of auc
  ci.list <- lapply(rocobj.list, function(x){ci.auc(x, conf.level = 0.95, method = "bootstrap")[c(1,3)]})
 
  # Comapre the AUC and calculate p value
  compare.p.value <- unlist(lapply(rocobj.list[-1], function(x){roc.test(rocobj.list[[1]], x, method = "delong")$p.value}))
  
  if(is.null(col)){
    colors <- rainbow(length(rocobj.list))
  } else {
    colors <- col
  }
   
  # plot roc curve
  if(is.null(pdf.file)){
    plot.roc(rocobj.list[[1]], legacy.axes = T, col = colors[1])
    for(i in 2:length(rocobj.list)){
      lines.roc(rocobj.list[[i]], legacy.axes = T, col = colors[i])
    }
    legend("bottomright", legend = names(rocobj.list), col = colors, lwd = 2)
  } else {
    pdf(pdf.file)
    plot.roc(rocobj.list[[1]], legacy.axes = T, col = colors[1])
    for(i in 2:length(rocobj.list)){
      lines.roc(rocobj.list[[i]], legacy.axes = T, col = colors[i])
    }
    legend("bottomright", legend = names(rocobj.list), col = colors, lwd = 2)
    dev.off()    
  }

  # return the main result
  model.names <- names(models)
  AUC <- round(unlist(auc.list), digits = 4)
  CI  <- unlist(lapply(ci.list, function(x){paste(round(x[1], digits = 4), "-", round(x[2], digits = 4), sep = "")}))
  pvalue <- c("-", round(compare.p.value, digits = 4))

  result <- data.frame(AUC, CI, pvalue)
  names(result) <- c("AUC", "95% CI", "p-value")

  return(result)
}





#-------------------------------------------------------------------------------------
#'@description 计算两个模型比较的NRI和IDI值
#'@param clinical.data：数据框，包含生存模型中所需的变量信息
#'@param time:数值型向量，患者对应的生存时间
#'@param event: 患者对应的生存状态，通常0=alive, 1=dead
#'@param models：列表，每一个元素是一个字符型向量，包含一个model中所有变量的名称
#'@param upper.time: 数值，感兴趣事件截止的时间点
#'@param npert：扰动重抽样的次数，默认是300，如果设为0，则不能计算置信区间和p值
#'@return 返回一个包含两个元素的list，第一个元素是IDI的信息（包括：IDI值，置信区间，p值），第二个元素是NRI的信息（包括：NRI值，置信区间，p值）
NRI.IDI <- function(clinical.data, time, event, models, upper.time, npert = 300){
  suppressPackageStartupMessages(require(survival))
  suppressPackageStartupMessages(require(survIDINRI))

  ##-------------------------------
  # CoxphObject：create a coxph object
  CoxphObject <- function(data, time, event, variables){
    if (length(variables) == 1) {
	  formula <- as.formula(paste("Surv(time, event)~", variables))
	} else {
	  formula <- as.formula(paste("Surv(time, event)~", paste(variables, collapse=" + ")))
	}      
    coxph.object <-coxph(formula, data = data)

    return(coxph.object)
  }
  
  coxph.list <- lapply(models, function(x){CoxphObject(data = clinical.data, time = time, event = event, variables = x)})

  covs.list <- lapply(coxph.list, function(x){as.matrix(model.matrix(object = x, data = clinical.data))})
  
  res.IDI.INF.list <- lapply(covs.list[-1], function(x){IDI.INF(indata = data.frame(time = time, event = event),
                                                        covs0 = covs.list[[1]],          
                                                        covs1 = x, 
                                                        t0 = upper.time,
                                                        npert = npert, npert.rand = NULL, seed1 = 1234)})

  ##------------------------------------------
  # 改造原来的IDI.INF.OUT，使它单独输出IDI或者NRI
  IDI.OUT <- function(x){
    if(!is.null(x$m1)){
      tmp <- x$m1
      names(tmp) <- c("IDI","Lower", "Upper","p-value")
      return(round(tmp, digits = 4))
    }else{
      tmp <- x$m1.est[1]
      names(tmp) <- "IDI"
      return(round(tmp, digits = 4))
    }
  }

  NRI.OUT <- function(x){
    if(!is.null(x$m2)){
      tmp <- x$m2
      names(tmp) <- c("NRI","Lower", "Upper","p-value")
      return(round(tmp, digits = 4))
    }else{
      tmp <- x$m2.est[1]
      names(tmp) <- "NRI"
      return(round(tmp, digits = 4))
    }
  }
  
  IDI.result.list <- lapply(res.IDI.INF.list, function(x){IDI.OUT(x)})
  NRI.result.list <- lapply(res.IDI.INF.list, function(x){NRI.OUT(x)})

  IDI.result <- do.call('rbind', IDI.result.list)
  NRI.result <- do.call('rbind', NRI.result.list)
  IDI.conf.interval <- paste(IDI.result[, 'Lower'], IDI.result[, 'Upper'], sep = "-")
  NRI.conf.interval <- paste(NRI.result[, 'Lower'], NRI.result[, 'Upper'], sep = "-")

  IDI.result <- data.frame(IDI.result[, 'IDI'], IDI.conf.interval, IDI.result[, 'p-value'])
  NRI.result <- data.frame(NRI.result[, 'NRI'], NRI.conf.interval, NRI.result[, 'p-value'])
  
  colnames(IDI.result) <- c("IDI", "95%CI", "p-value")
  colnames(NRI.result) <- c("NRI", "95%CI", "p-value")
  rownames(IDI.result) <- paste(rownames(IDI.result), " VS ", names(models)[1],sep = "")
  rownames(NRI.result) <- paste(rownames(NRI.result), " VS ", names(models)[1],sep = "")

  result <- list(IDI = IDI.result, NRI = NRI.result)
  return(result)
}





#----------------------------------------------------------------------------
#'@description 计算模型的AIC值
#'@param clinical.data：数据框，包含生存模型中所需的变量信息
#'@param event: 患者对应的生存状态，通常0=alive, 1=dead
#'@param time:数值型向量，患者对应的生存时间
#'@param models：列表，每一个元素是一个字符型向量，包含一个model中所有变量的名称
#'@return 返回一个矩阵，第一列代表模型的自由度（df），第二列代表模型的AIC值，每行对应一个模型的信息
aic <- function(clinical.data, time, event, models){
  suppressPackageStartupMessages(require(survival))

  ##-------------------------------
  # CoxphObject：create a coxph object
  CoxphObject <- function(data, time, event, variables){
    if (length(variables) == 1) {
	  formula <- as.formula(paste("Surv(time, event)~", variables))
	} else {
	  formula <- as.formula(paste("Surv(time, event)~", paste(variables, collapse=" + ")))
	}      
    coxph.object <-coxph(formula, data = data)

    return(coxph.object)
  }

  coxph.list <- lapply(models, function(x){CoxphObject(data = clinical.data, time = time, event = event, variables = x)})
  aic.list <- lapply(coxph.list, function(x){extractAIC(x)})
  aic.result <- do.call('rbind', aic.list)
  colnames(aic.result) <- c("df", "AIC")
    
  return(aic.result)
}




