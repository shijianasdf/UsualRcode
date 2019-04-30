# 两组数据的差异
# student's t.test判断两组数据均值差异是否显著   alternative:表示备择假设 x>y 默认two.sided
round(t.test(x,y,alternative="greater")$p.value,digits=4)
# wilcox.test判断两组数据总体分布是否存在差异  alternative:备择假设 x>y 默认two.sided
wilcox.test(x,y,alternative="greater")$p.value   #wilcox要求两组向量较长一些 

# 独立性检验
# fisher.exact 原假设:OR=1 备择假设；OR>1
a <- unclass(table(filtered.new.table[,i], filtered.new.table[,j]));#table接收2个因子,生成检验的表格
fisher.test(a)$p.value;
fisher.test(matrix(c(12,14,23,34),nrow=2),alternative="greater");
fisher.test(matrix(c(12,14,23,34),nrow=2),alternative="greater");
# chisq.test
chisq.test(matrix(c(12,14,23,34),nrow=2));
chisq.test(a)$p.value; #上面的表格

#相关性
cor(matrixdata);
cor.test(x,y); #x,y长度相同的两个向量
cor.test(~col1+col2,data=matrixdata);

# 累积超几何检验  m代表交集个数,M代表先验基因集合基因个数(go term等),N代表背景基因集合基因个数,n代表用户提供的基因集合基因个数
phyper(m-1,M,N-M,n,lower.tail=F,log.p=F); #看交集是否显著
phyper(m,M,N-M,n,lower.tail=T,log.p=F); #看互斥是否显著
# 累积二项分布检验  m次伯努利二项实验(只有两个结果的实验，每次成功的概率为p，例如抛硬币)，
pbinom(n-1,m,p,lower.tail=F);
# 二项分布检验
binom.test(n,m,p);
# 累积泊松分布检验 染色体上一段固定区域的突变个数服从泊松分布，如果落在长度为10bp区域的突变个数为n个，那么n个突变是否显著的大
ppois(n-1,lamada,lower.tail=F);
# lamada <- 10bp区域内突变个数的均值

#生存分析
library(splines);
library(survival);
# Create the simplest test data set
test1 <- list(time=c(4,3,1,1,2,2,3),

              status=c(1,1,1,0,1,1,0),

              x=c(0,2,1,1,1,0,0),

              sex=c(0,0,0,0,1,1,1));

# Fit a stratified model
coxph(Surv(time, status) ~ x + strata(sex), test1); #strata函数对性别的这个混杂因素进行分层，#我理解就是分组，考虑在sex分组情况下x因素对生存长短的影响。coxph(Surv(time, #status) ~ x + sex, test1)，这个是分别考虑x，sex因素对生存长短的影响。
head(data);
#   ID Clinic Status Days.survival Prison Dose
# 1  1      1      1           428      0   50
# 2  2      1      1           275      1   55
# 3  3      1      1           262      0   55
# 4  4      1      1           183      0   30
# 5  5      1      1           259      1   65
# 6  6      1      1           714      0   55
#构建生存对象   Surv(生存时间，生存状态)
surv <- Surv(data$Days.survival,data$Status);
surv
 #  [1]  428   275   262   183   259   714   438   796+  892   393   161+  836
 # [13]  523   612   212   399   771   514   512   624   209   341   299   826+
 # [25]  262   566+  368   302   602+  652   293   564+  394   755   591   787+
kmfit <- survfit(surv~data$Prison); #根据prison分组绘制K-M生存曲线
p_value <- pchisq(survdiff(y~t(tempMatrix2[i,]))$chisq, 1, lower.tail = FALSE); #log rank检验
plot(kmfit1,lty = c('solid', 'dashed'),col=c('red','blue'),xlab='survival time in days',ylab='survival probabilities',main=p_value,sub=names[i]);
legend('topright', c('Clinic 1','Clinic 2'), lty=c('solid','dashed'), col=c('red','blue'));

#拟合优度检验
#Kolmogorov-Smirnov是比较一个频率分布f(x)与理论分布g(x)或者两个观测值分布的检验方法。
#其原假设H0:两个数据分布一致或者数据符合理论分布。
#D=max| f(x)- g(x)|，当实际观测值D>D(n,α)则拒绝H0，否则则接受H0假设。
#R语言中ks.test有四个参数，第一个参数x为观测值向量，
#第二个参数y为第二观测值向量或者累计分布函数或者一个真正的累积分布函数如pnorm，只对连续CDF有效。
#第三个参数为指明是单侧检验还是双侧检验，exact参数为NULL或者一个逻辑值，表明是否需要计算精确的P值。
ks.test(rnorm(100),rnorm(50)) #两个观测向量
ks.test(rnorm(100),"pnorm") #一个观测向量和一个理论分布
ks.test(rnorm(100),"punif") #一个观测向量和一个理论分布


