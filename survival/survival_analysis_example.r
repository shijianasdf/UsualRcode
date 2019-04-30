# 首先需要构造生存分析表
# 以gene.sur.data.clonality为例
# patient Clonal_Status OS_MONTHS OS_STATUS(数值型) AGE tumor_purity race ethnicity IDH_STATUS
#  8864        1            12         1     45      0.7      **      **        1

# 导入生存分析包
library(survival);
##生存对象构建
y <- with(gene.sur.data.clonality, Surv(OS_MONTHS, OS_STATUS));
## log-rank检验
clonality.dfs.result <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ Clonal_Status, data = gene.sur.data.clonality); # log-rank检验
clonality.dfs.p.value <- 1 - pchisq(clonality.dfs.result$chisq, length(clonality.dfs.result$n) - 1); # 

#计算log-rank检验p值
library(survminer)
#生成log-rank检验的p值,保留两位有效数字
pval <- round(surv_pvalue(km.curves, gene.sur.data.clonality)$pval,digits=4);
#计算pairwise的log rank的p值
res <- pairwise_survdiff(Surv(time, event)~sample.label, data=gene.sur.data.clonality);
#计算中位生存期以及CI
median.os <- surv_median(fit=km.curves);

## 生存曲线绘图  根据Clonal_status给样本贴标签,Clonal_status可以是任意值(数字，字符串)
km.curves <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ Clonal_Status, data=gene.sur.data.clonality);   
plot(km.curves, lty = c('solid', 'dashed'), col= c('red','blue'), xlab = 'survival time in months', ylab = 'survival probabilities', main = gene.name, bty = 'n');
legend('topright', c('Clonal','Subclonal'), lty = c('solid','dashed'), col = c('red','blue'), bty = 'n');
legend('bottomleft', paste("p=", as.character(round(clonality.dfs.p.value, 6))), bty = 'n');

library(survminer)
ggsurv <- ggsurvplot(                         #注意此时并没有画图，只是存在ggsurv变量中
                         km.curves,               # survfit object with calculated statistics.
                         data = gene.sur.data.clonality,             # data used to fit survival curves.
                         palette = c("red","blue"), #分组颜色
                         
                         #图的主题构架
                         risk.table = T,       # 是否展示风险table.
                         pval = T,             # 是否展示log-rank检验的p值.
                         surv.median.line = c("none", "hv", "h", "v"),  # 是否展示中位生存期.surv_median
						 conf.int = T, #是否画出置信区间
                         title = "shijian",     #主标题名字
                         font.main = 15,       #主标题字体大小              
                         xlab = paste("Time","in","Monthes",sep = " "),   # customize X axis label.
                         ylab = "Overall Survival",   # customize Y axis label.
                         
                         #图例设置
                         legend.title = "", #图例标题，一般不用，设为空
                         legend.labs = substr(names(km.curves$strata),start = 14,stop = 1000), #图例文字描述
                         legend = c(0.8,0.9), #图例的位置，取值在【0,1】之间
                         font.legend = 9,     #图例字体大小
                         
                         #risk table设置
                         tables.theme = theme_cleantable(),#table主题
                         risk.table.title = "No. at risk:",#table标题
                         risk.table.y.text.col = T, # 使用颜色代替Y轴文字
                         risk.table.y.text = FALSE, # Y轴不使用文字注释
                         tables.height = 0.15,      # table的高度
                         risk.table.fontsize = 3,    # risk table内文字的大小
						             ggtheme = theme_bw() # Change ggplot2 theme
                        );
    ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"));#KM曲线图的标题居中
    ggsurv$table <- ggsurv$table + theme(plot.title = element_text(hjust = -0.04), plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"));#风险table图的标题居中
	
## 计算1,2,3,4,5年生存率
survival_rate <- summary(km.curves,times=c(1,2,3,4,5))

# Cox回归检验
# 单因素回归检验
coxmodel <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ Clonal_Status, data = gene.sur.data.clonality);
summary(coxmodel); #coxph默认将离散型变量转换为因子，并将level中的第一个作为HR参考项.
# 多因素回归检验 #默认将离散型变量转换为因子，并将level中的第一个作为HR参考项.
coxmodel <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ Clonal_Status + AGE + tumor_purity + race + ethnicity + IDH_STATUS, data = gene.sur.data.clonality); #默认将离散型变量转换为因子，并将level中的第一个作为HR参考项.
summary(coxmodel);
# 互作COX模型
coxmodel <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ Clonal_Status + AGE + tumor_purity + race + ethnicity*IDH_STATUS, data = gene.sur.data.clonality);
summary(coxmodel);

#cox回归模型解释：
#1、统计学显著性
#z值是由wald统计量，z = coef/se(coef) ， 而wald统计检验估计的是 是否一个给定变量的系数 beta 不等于0是显著的。由以上结果可以看出P=0.00149，说明sex具有统计学显著的系数。
#2、回归系数coefs
#如果回归系数为正值，说明死亡风险高，预后差；为负值，则风险低，预后好。
#变量sex为数值型变量1ale 2emale，HR则是第二组相对于第一组来说的，也就是female versus male. 在此beta系数为-0.531说明female比male有更低的死亡风险(低生存率)。
#3、风险比HR=exp(coef)
#风险比，给定了sex的影响大小。例如，female 降低了41%的风险，female与好的预后相关
#For example, being female (sex=2) reduces the hazard by a factor of 0.59, or 41%. Being female is associated with good prognostic.
#不过我们通常会对连续变量进行单因素cox回归分析，则只需要关注HR是否大于1，大于1为风险因素，小于1为保护因素。
#4、HR的置信区间
#结果中给出的HR(exp(coefs))的置信区间：lower 95% bound = 0.4237, upper 95% bound = 0.816
#5、该模型的整体的统计学显著性
#这个output中给出了三个检验模型的P值：The likelihood-ratio test, Wald test, and score logrank #statistics.这三种方法是类似的，如果样本量足够大，则三个方法会得到相同的结果，当样本量较小时，The #likelihood-ratio test能得到较好的结果。 


res.cox <- step(coxph(Surv(time, status) ~ age + sex + ph.ecog, data =  lung),direction = "backward")
#####################利用maxstat判断分组阈值，以达到最好预后生存预测##############################
library("maxstat");
library("survival");
data("DLBCL", package="maxstat");
head(DLBCL);
# RowNames    DLCLid              GEG time cens IPI          MGE
# 1        7 DLCL-0001        GC B-Like 77.4    0   2  0.247973684
# 2        8 DLCL-0002 Activated B-like  3.4    1   3 -0.006894737
# 3        9 DLCL-0003        GC B-Like 71.3    1   2  0.137552632
# 4       10 DLCL-0004        GC B-Like 69.6    0   0  0.371289474
# 5       11 DLCL-0005 Activated B-like 51.2    0   1  0.025710526
# 6       12 DLCL-0006 Activated B-like  3.2    1   3  0.098131579
mtHL <- maxstat.test(Surv(time, cens) ~ MGE,
        data=DLBCL, smethod="LogRank", pmethod="HL");
mtHL
#Maximally selected LogRank statistics using HL
#data:  Surv(time, cens) by MGE
#M = 3.171, p-value = 0.02218
#sample estimates:
#estimated cutpoint 
#         0.1860526
