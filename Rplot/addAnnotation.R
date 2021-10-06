library(ggplot2)
library(ggannotate)

p <- ggplot(mtcars, 
            aes(x = wt, y = mpg)) + 
  geom_point() 

ggannotate(p)
