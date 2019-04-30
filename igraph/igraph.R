library(igraph)
head(e.p)
  enhancer       gene.p
1     Enh1 RP11-34A14.3
2     Enh1        LOXL4
3     Enh1         HPS1
4     Enh1        HPSE2
5     Enh1      ZDHHC16
6     Enh1   AL355490.2
e.p <- na.omit(e.p)
g <- graph_from_data_frame(e.p, directed=FALSE) #data.frame to igraph 构建igraph对象,directed 参数控制graph 有无方向
#V(g)和E(g)可以用来查看网络g的节点和边  
V(g)  V(g)$name #点的名字
E(g)  E(g)[node1%--%node2]  #node1对应node2的边 %--%代表无向网络
#提取最大成分子网
max.g <- clusters(g)
sub.g <- components(g)
subnetList <- groups(sub.g) #各个子网中的点 
max.son <- subnetList[[which.max(sapply(subnetList,function(x){length(x)}))]]
E(g)[max.son%--%max.son] #取边
max.subnet <- induced_subgraph(g,max.son) #从g中提取max.son子节点构建子网igraph对象
plot(max.subnet,layout=layout_with_kk, vertex.color="green")
max.subnet <- as.data.frame(get.edgelist(max.subnet)) #get.edgelist(max.subnet)将igraph对象转换为矩阵,再转换为data.frame
max.subnet <- graph_from_data_frame(max.subnet,directed=FALSE) #data.frame转换为igraph对象
degree(max.subnet) #求该网的度
closeness(max.subnet,V(g))




