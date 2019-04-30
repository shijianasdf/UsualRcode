library(randomForestSRC)
data(pbc, package = "randomForestSRC")
pbc.obj2 <- rfsrc(Surv(days, status) ~ ., pbc,
                  nsplit = 10, na.action = "na.impute")
print(pbc.obj2)
plot(pbc.obj2)
