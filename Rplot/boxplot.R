###################################
#boxplot单组或者多组数据的分布图
###################################
boxplot(x, ..., range = 1.5, width = NULL, varwidth = FALSE,
        notch = FALSE, outline = TRUE, names, plot = TRUE,
        border = par("fg"), col = NULL, log = "",
        pars = list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5),
        horizontal = FALSE, add = FALSE, at = NULL)
#x: 列表或者向量
#outline: 是否画离群点

######################################
#
#######################################
ggplot(data=mpg,mapping=aes(x=manufacturer,y=cty))+
geom_boxplot(mapping=aes(colour=manufacturer))+
theme_bw()