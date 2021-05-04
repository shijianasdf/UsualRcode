#重启动随机游走
install.packages("RANKS")
install.packages("bionetdata")
library(RANKS)
library(bionetdata);
data(DD.chem.data);
data(DrugBank.Cat);
labels <- DrugBank.Cat[,"Penicillins"];
ind.pos <- which(labels==1);
# 2-step RWR
res <- RWR(DD.chem.data, ind.pos, tmax = 2);
## Not run: 
# till to convergence
res <- RWR(DD.chem.data, ind.pos, tmax = 5000, eps=1e-6);
# 5 steps and higher gamma
res <- RWR(DD.chem.data, ind.pos, tmax = 5, gamma=0.8);

#创建网络,根据网络获取邻接矩阵
library(igraph)
data_stru<-graph.data.frame(data)
get.adjacency(data_stru, attr="weight")
