data <- read.table('X:\\bipartite_network.txt',header = T)

#创建网络,画图
library(igraph)
data_stru<-graph.data.frame(data)
plot(data_stru,edge.width=2,edge.arrow.size=0.5,edge.color=c(data$OUT))

#网络的拓扑属性
summary(data_stru) #简单查看网络节点数、边数等
V(data_stru)  #查看网络节点
V(data_stru)[degree(data_stru)>4] #查看节点数大于4的节点
E(data_stru)  #查看网络的边
degree(data_stru) #网络中每个节点度大小
degree.distribution(data_stru)  #计算节点度分布（度=c(0,1,2,3...)的节点数占总节点的比例）
transitivity(data_stru) #计算整个网络的聚集系数
transitivity(data_stru,"local") #每个节点的聚集系数
transitivity(data_stru,"local") #计算每个节点的聚集系数
diameter(data_stru) #计算网路的直径，即MAX（任意两个节点之间的最短路径）
average.path.length(data_stru)  #计算整个网路的平均距离（任意两节点直接的最短路径的平均值）
clusters(data_stru)$no  #查看孤立点个数
closeness(data_stru)  #计算紧密度
betweenness(data_stru)  # 计算介质中心性
is.directed(data_stru)  #网络是否为有向图
is.connected(data_stru) #网络是否连通
no.clusters(data_stru) #网络有多少分支
graph.density(data_stru)  #网络的密度
#以下函数中‘in’,‘out’为有向图的概念，无向图去掉mode参数
max(degree(data_stru, mode="in")) #有向图中in节点度最大为多少
max(degree(data_stru, mode="out"))  #有向图中out节点度最大为多少
plot(degree.distribution(data_stru, mode="in")) #画in节点度分布图
plot(degree.distribution(data_stru, mode="in"), log="xy") #画对数横纵坐标的in节点度分布图
plot(degree.distribution(data_stru, mode="out"), log="xy")  #画对数横纵坐标的out节点度分布图

# 构建一个与原网络有相同度数和相同节点数的随机网络 
#方法一
g <- erdos.renyi.game(vcount(data_stru), ecount(data_stru), type="gnm") 
plot(g,edge.width=2,edge.arrow.size=0.5)
#方法二
g2 <- degree.sequence.game(degree(data_stru,mode="all"), method="vl") 
plot(g2,edge.width=2,edge.arrow.size=0.5)

