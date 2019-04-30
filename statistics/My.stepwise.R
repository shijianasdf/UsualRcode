library(My.stepwise)
library(help="My.stepwise")
#usage:My.stepwise.glm(Y, variable.list, in.variable = "NULL", data, sle = 0.15,
#                sls = 0.15, myfamily, myoffset = "NULL")
my.variable.list=c("SEX","MSS_state","AJCC_TUMOR_PATHOLOGIC_PT","AJCC_NODES_PATHOLOGIC_PN","AJCC_PATHOLOGIC_TUMOR_STAGE","RESIDUAL_TUMOR","AGE","subtype","tumour_site_2")
My.stepwise.glm(Y="AJCC_METASTASIS_PATHOLOGIC_PM",
                variable.list = my.variable.list,
                in.variable=c("APC"),
                myfamily="binomial",
                data=logist.clinical.M.data[[1]])
head(logist.clinical.M.data[[1]])
colnames(logist.clinical.M.data[[1]])

data("iris")
names(iris)
my.data <- iris[51:150, ]
my.data$Width <- (my.data$Sepal.Width + my.data$Petal.Width)/2
names(my.data)
dim(my.data)
my.data$Species1 <- ifelse(my.data$Species == "virginica", 1, 0)
my.variable.list <- c("Sepal.Length", "Petal.Length")
My.stepwise.glm(Y = "Species1", variable.list = my.variable.list,
                in.variable = c("Width"), data = my.data, myfamily = "binomial")

my.variable.list <- c("Sepal.Length", "Sepal.Width", "Width")
My.stepwise.glm(Y = "Species1", variable.list = my.variable.list,
                data = my.data, sle = 0.25, sls = 0.25, myfamily = "binomial")
