##提取COX模型的HR,CI,p值
  mvcoxres <- function(mod){
    hr.ci <- paste0(
      round(summary(mod)$conf.int[, 1], 2), " (",
      round(summary(mod)$conf.int[, 3], 2), ", ",
      round(summary(mod)$conf.int[, 4], 2), ")"
    )
    p <- round(summary(mod)$coefficients[, 5], 3)
    res <- data.frame(hr.ci, p, stringsAsFactors = F)
    res$p[res$p < 0.001] <- "<.001"
    colnames(res) <- c("**HR (95% CI)**", "**p-value**")
    return(res)
  }
  library(MASS)
  library(survival)
  #pos <- complete.cases(COX.OS.all.gene[,c(2,3)])
  head(COX.OS.all.gene)
  coxmodel <- coxph(as.formula(paste("Surv(time,event)~",paste(colnames(COX.OS.all.gene)[-c(1:3)],collapse = "+"))),
                    data=na.omit(COX.OS.all.gene)) #stepAIC 不允许生存数据有NA
  summary(coxmodel)
  extractAIC(coxmodel)
  step.cox <- stepAIC(coxmodel,direction = "both",trace = F) #stepAIC 不允许生存数据有NA
  #step.cox.result <- summary(step.cox)
  stepwise.table <- mvcoxres(step.cox)
