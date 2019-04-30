####--------------------------------------
#Taken together, the main purpose of principal component analysis is to:
#   identify hidden pattern in a data set,
#   reduce the dimensionnality of the data by removing the noise and redundancy in the data,
#   identify correlated variables
####--------------------------------------
library("FactoMineR")
library("factoextra")
library(pca3d)
data(decathlon2)
head(decathlon2)
decathlon2.active <- decathlon2[1:23, 1:10]
head(decathlon2.active[, 1:6], 4)
res.pca <- PCA(decathlon2.active, graph = FALSE)
print(res.pca)
eig.val <- get_eigenvalue(res.pca)
eig.val
var <- get_pca_var(res.pca)
var
# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)
fviz_pca_var(res.pca, col.var = "black")

