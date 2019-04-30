library(igraph)

networkPlot <- function(x){   #x类型data.frame,网络关系对数据，含两列，分别是节点1和节点2。
  g <- graph.data.frame(x)    #生成网络。
  
  node1 <- table(degree(g)[1:length(levels(factor(x[[1]])))])    #查看节点1度分布。
  Node1 <- rbind(as.numeric(names(node1)),node1)                 #把节点1分布数值化。
  m_lm1 <- lm(log10(Node1[2,])~log10(Node1[1,]))                 #度和节点数都取对数后做线性回归。
  
  node2 <- table(degree(g)[-(1:length(levels(factor(x[[1]]))))]) #节点2处理方法同节点1。
  Node2 <- rbind(as.numeric(names(node2)),node2)
  m_lm2 <- lm(log10(Node2[2,])~log10(Node2[1,]))
  pdf("a:/DegreeDistribution.pdf")  #画PDF图。
  #--------------------------------------------------------------------------------------------------
  #节点1度分布图（x中的第一列）
  plot(Node1[1,],Node1[2,],las=1,xlab='Degree',ylab='Number of Node1')    #画度及其节点数的散点图
  lines(Node1[1,],10^( m_lm1$coefficients[1])*Node1[1,]^(m_lm1$coefficients[2]),lwd=2)   #添加回归函数线。
  #节点1度分布图，使用对数横纵坐标
  plot(Node1[1,],Node1[2,],log='xy',las=1,xlab='Degree',ylab='Number of Node1')
  lines(Node1[1,],10^( m_lm1$coefficients[1])*Node1[1,]^(m_lm1$coefficients[2]),lwd=2)
  #--------------------------------------------------------------------------------------------------
  #节点2度分布（x中的第二列
  plot(Node2[1,],Node2[2,],las=1,xlab='Degree',ylab='Number of Node2')
  lines(Node2[1,],10^( m_lm2$coefficients[1])*Node2[1,]^(m_lm2$coefficients[2]),lwd=2)
  #节点2度分布（文件中的第二列
  plot(Node2[1,],Node2[2,],log='xy',las=1,xlab='Degree',ylab='Number of Node2')
  lines(Node2[1,],10^( m_lm2$coefficients[1])*Node2[1,]^(m_lm2$coefficients[2]),lwd=2)
  #--------------------------------------------------------------------------------------------------
  dev.off()
}

#networkPlot(targets)
