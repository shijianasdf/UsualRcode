####################################################ROC###################################################
library(pROC);
#Tools for visualizing, smoothing and comparing receiver operating characteristic (ROC curves).
#(Partial) area under the curve (AUC) can be compared with statistical tests based on U-statistics or
#bootstrap. Confidence intervals can be computed for (p)AUC or ROC curves. Sample size / power
#computation for one or two ROC curves are available.
#真实类别(1代表case,0代表control)
obs<-rep(0:1, each=50); 
#预测属于1case的概率
pred<-c(runif(50,min=0,max=0.8),runif(50,min=0.3,max=0.6)); 
pred1<-c(runif(50,min=0,max=0.8),runif(50,min=0.3,max=0.6));
#roc对象
roc <- roc(obs,pred);
#smooth.roc对象
smooth.roc <- roc(obs,pred1,smooth=T); smooth.roc <- smooth(roc);
#绘制ROC曲线,还有很多有用参数，具体自己查
pdf();
plot.roc(smooth.roc,xlim=c(0,1),ylim=c(0,1)); 
lines(roc1); #在上一幅图中再次加入ROC曲线
dev.off();
#计算AUC
auc <- auc(roc); 
#通过bootstrap方式计算得到CI
CI <- ci(roc); 
#通过bootstrap方式计算得到CI.auc
CI.auc <- ci.auc(roc); 
CI.se <- ci.se(roc); 
#x:best all local maximas返回ROC上面的坐标点(specificity,sensitivity)，其中best返回最优cutoff，方法为Youden index
coords <- coords(roc,x="best"); 