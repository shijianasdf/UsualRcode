require(ggplot2)
data(diamonds)
set.seed(42)
small <- diamonds[sample(nrow(diamonds), 1000), ]
head(small)

theme_set(theme_bw()); #白色背景
p1 <- ggplot(data=small)+geom_point(mapping=aes(x=carat,y=price,shape=cut, colour=color));
p2 <- ggplot(data=small)+geom_boxplot(aes(x=cut, y=price,fill=color));
require(ggpubr); #布局
ggarrange(p1, p2,nrow=2,heights=c(1,1));