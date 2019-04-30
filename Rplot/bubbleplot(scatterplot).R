rm(list = ls())
library(datasets)
library(ggplot2)
data(airquality)
aq_trim <- airquality[which(airquality$Month == 7 |
                              airquality$Month == 8 |
                              airquality$Month == 9), ]
aq_trim$Month <- factor(aq_trim$Month,
                        labels = c("July", "August", "September"))
fill = c("steelblue", "yellowgreen", "violetred1")
# Ozone Solar.R Wind Temp     Month Day
# 62    135     269  4.1   84      July   1
# 63     49     248  9.2   85      July   2
# 64     32     236  9.2   81      July   3
# 65     NA     101 10.9   84      July   4
# 66     64     175  4.6   83      July   5

p6 <- ggplot(aq_trim, aes(x = Day, y = Ozone, size = Wind, fill = Month)) +
  geom_point(shape = 21) +
  ggtitle("Air Quality in New York by Day") +
  labs(x = "Day of the month", y = "Ozone (ppb)",
       size = "Wind Speed (mph)", fill = "Months") +
  scale_x_continuous(breaks = seq(1, 31, 5)) +
  scale_size(range = c(1, 10)) +
  scale_fill_manual(values = fill) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.size = unit(1, "cm"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 9),
        axis.text.y=element_text(colour="black", size = 9))
p6
