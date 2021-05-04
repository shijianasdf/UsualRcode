#重启动随机游走
install.packages("RANKS")
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
