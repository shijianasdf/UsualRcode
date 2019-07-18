#----------------------
#' combine plot
#----------------------

#方法1  par(mfrow=c(2,2)) par(mfcol=c(2,2)) 
{
  #Graphical parameter mfrow can be used to specify the number of subplot we need.
  #It takes in a vector of form c(m, n) which divides the given plot into m*n array of subplots. 
  #For example, if we need to plot two graphs side by side, we would have m=1 and n=2. Following example illustrates this.
  opar <- par(no.readonly = T)
  par(mfrow=c(1,2))    # set the plotting area into a 1*2 array
  max.temp <- c(22,27,26,24,23,26,28)
  names(max.temp) <- c("Sun", "Mon", "Tue", "Wen", "Thu", "Fri", "Sat")
  barplot(max.temp, main="Barplot")
  pie(max.temp, main="Piechart", radius=1)
  par(opar)
  
  Temperature <- airquality$Temp
  Ozone <- airquality$Ozone
  par(mfrow=c(2,2))
  hist(Temperature)
  boxplot(Temperature, horizontal=TRUE)
  hist(Ozone)
  boxplot(Ozone, horizontal=TRUE)
  par(opar)
  
  #Same plot with the change par(mfcol = c(2, 2)) would look as follows. 
  #Note that only the ordering of the subplot is different.
  Temperature <- airquality$Temp
  Ozone <- airquality$Ozone
  par(mfcol = c(2, 2))
  hist(Temperature)
  boxplot(Temperature, horizontal=TRUE)
  hist(Ozone)
  boxplot(Ozone, horizontal=TRUE)
  par(opar)
  
  #The graphical parameter fig lets us control the location of a figure precisely in a plot.
  #We need to provide the coordinates in a normalized form as c(x1, x2, y1, y2). 
  #For example, the whole plot area would be c(0, 1, 0, 1) with (x1, y1) = (0, 0) being the lower-left corner and (x2, y2) = (1, 1) being the upper-right corner.
  
  # make labels and margins smaller
  par(cex=0.7, mai=c(0.1,0.1,0.2,0.1))
  Temperature <- airquality$Temp
  # define area for the histogram fig= starts a new plot, so to add to an existing plot use new=TRUE.
  par(fig=c(0.1,0.7,0.3,0.9))
  hist(Temperature)
  # define area for the boxplot
  par(fig=c(0.8,1,0,1), new=TRUE)
  boxplot(Temperature)
  # define area for the stripchart
  par(fig=c(0.1,0.67,0.1,0.25), new=TRUE)
  stripchart(Temperature, method="jitter")
  par(opar)
}

#方法2 layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE)) 
{
  # One figure in row 1 and two figures in row 2
  attach(mtcars)
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  hist(wt)
  hist(mpg)
  hist(disp)
  
  # One figure in row 1 and two figures in row 2
  # row 1 is 1/3 the height of row 2
  # column 2 is 1/4 the width of the column 1 
  attach(mtcars)
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
         widths=c(3,1), heights=c(1,2))
  hist(wt)
  hist(mpg)
  hist(disp)
  
  # Add boxplots to a scatterplot fig= starts a new plot, so to add to an existing plot use new=TRUE.
  par(fig=c(0,0.8,0,0.8), new=TRUE)
  plot(mtcars$wt, mtcars$mpg, xlab="Car Weight",
       ylab="Miles Per Gallon")
  par(fig=c(0,0.8,0.55,1), new=TRUE)
  boxplot(mtcars$wt, horizontal=TRUE, axes=FALSE)
  par(fig=c(0.65,1,0,0.8),new=TRUE)
  boxplot(mtcars$mpg, axes=FALSE)
  mtext("Enhanced Scatterplot", side=3, outer=TRUE, line=-3)
 
}

#方法3  gridExtra grid
#To arrange multiple ggplot2 graphs on the same page, 
#the standard R functions - par() and layout() - cannot be used.
#The basic solution is to use the gridExtra R package, which comes with the following functions:
#grid.arrange() and arrangeGrob() to arrange multiple ggplots on one page
#marrangeGrob() for arranging multiple ggplots over multiple pages.
{
  library(ggpubr)
  # ToothGrowth
  data("ToothGrowth")
  head(ToothGrowth)
  # mtcars 
  data("mtcars")
  mtcars$name <- rownames(mtcars)
  mtcars$cyl <- as.factor(mtcars$cyl)
  head(mtcars[, c("name", "wt", "mpg", "cyl")])
  # Box plot (bp)
  bxp <- ggboxplot(ToothGrowth, x = "dose", y = "len",
                   color = "dose", palette = "jco")
  bxp
  # Dot plot (dp)
  dp <- ggdotplot(ToothGrowth, x = "dose", y = "len",
                  color = "dose", palette = "jco", binwidth = 1)
  dp
  # Bar plot (bp)
  bp <- ggbarplot(mtcars, x = "name", y = "mpg",
                  fill = "cyl",               # change fill color by cyl
                  color = "white",            # Set bar border colors to white
                  palette = "jco",            # jco journal color palett. see ?ggpar
                  sort.val = "asc",           # Sort the value in ascending order
                  sort.by.groups = TRUE,      # Sort inside each group
                  x.text.angle = 90           # Rotate vertically x axis texts
  )
  bp + font("x.text", size = 8)
  # Scatter plots (sp)
  sp <- ggscatter(mtcars, x = "wt", y = "mpg",
                  add = "reg.line",               # Add regression line
                  conf.int = TRUE,                # Add confidence interval
                  color = "cyl", palette = "jco", # Color by groups "cyl"
                  shape = "cyl"                   # Change point shape by groups "cyl"
  )+
    stat_cor(aes(color = cyl), label.x = 3)       # Add correlation coefficient
  sp
  library("gridExtra")
  grid.arrange(bxp, dp, bp + rremove("x.text"), 
               ncol = 2, nrow = 2)
  
  grid.arrange(sp,                             # First row with one plot spaning over 2 columns
               arrangeGrob(bxp, dp, ncol = 2), # Second row with 2 plots in 2 different columns
               nrow = 2)                       # Number of rows
  
  grid.arrange(bp,                                    # bar plot spaning two columns
               bxp, sp,                               # box plot and scatter plot
               ncol = 2, nrow = 2, 
               layout_matrix = rbind(c(1,1), c(2,3)))
  
  library("gridExtra")
  library("cowplot")
  # Arrange plots using arrangeGrob
  # returns a gtable (gt)
  gt <- arrangeGrob(bp,                               # bar plot spaning two columns
                    bxp, sp,                               # box plot and scatter plot
                    ncol = 2, nrow = 2, 
                    layout_matrix = rbind(c(1,1), c(2,3)))
  # Add labels to the arranged plots
  p <- as_ggplot(gt) +                                # transform to a ggplot
    draw_plot_label(label = c("A", "B", "C"), size = 15,
                    x = c(0, 0, 0.5), y = c(1, 0.5, 0.5)) # Add labels
  p
  
  
  # The grid R package can be used to create a complex layout with the help of the function grid.layout(). It provides also the helper function viewport() to define a region or a viewport on the layout. The function print() is used to place plots in a specified region.
  # 
  # The different steps can be summarized as follow :
  #   
  # Create plots : p1, p2, p3, ….
  # Move to a new page on a grid device using the function grid.newpage()
  # Create a layout 2X2 - number of columns = 2; number of rows = 2
  # Define a grid viewport : a rectangular region on a graphics device
  # Print a plot into the viewport
  library(grid)
  # Move to a new page
  grid.newpage()
  # Create layout : nrow = 3, ncol = 2
  pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 2)))
  # A helper function to define a region on the layout
  define_region <- function(row, col){
    viewport(layout.pos.row = row, layout.pos.col = col)
  } 
  # Arrange the plots
  print(sp, vp = define_region(row = 1, col = 1:2))   # Span over two columns
  print(bxp, vp = define_region(row = 2, col = 1))
  print(dp, vp = define_region(row = 2, col = 2))
  print(bp + rremove("x.text"), vp = define_region(row = 3, col = 1:2))
}

# gridExtra详解
{
  library(gridExtra)
  library(ggplot2)
  library(lattice)
  
  p <- qplot(1,1)
  p2 <- xyplot(1~1)   ##lattice包
  
  grid.arrange(p,p2,ncol = 2)
  
  
  library(gridExtra)
  library(ggplot2)
  library(lattice)
  
  gs <- list(NULL)
  gs[[1]] <- qplot(1,1)
  gs[[2]] <- xyplot(1~1)   ##lattice包
  
  grid.arrange(grobs=gs,ncol = 2)
  
  lay <- rbind(c(1,1,1,2,3),
               c(1,1,1,4,5),
               c(6,7,8,9,9))
  
  grid.arrange(grobs = gs,layout_matrix = lay)
  
  library(gridExtra)
  library(ggplot2)
  
  g <- ggplotGrob(qplot(1, 1) +
                    theme(plot.background = element_rect(colour = "black")))
  qplot(1:10, 1:10) +
    annotation_custom(
      grob = g,
      xmin = 1, xmax = 5, ymin = 5, ymax = 10
    ) 
}


# DEMO
{
  library(ggplot2)
  library(gridExtra)
  ToothGrowth$dose <- as.factor(ToothGrowth$dose)
  head(ToothGrowth)
  g1 <- ggplot(ToothGrowth, aes(x=dose, y=len)) + 
          geom_boxplot()+stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
          geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
  g2 <- ggplot(ToothGrowth, aes(x=dose, y=len)) + 
          geom_violin()+stat_summary(fun.y=mean, geom="point", shape=23, size=2)+ 
          geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
  g3 <- ggplot(ToothGrowth, aes(x=dose, y=len)) + 
          geom_boxplot()+stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
          geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
  g4 <- ggplot(ToothGrowth, aes(x=dose, y=len)) + 
          geom_violin()+stat_summary(fun.y=mean, geom="point", shape=23, size=2)+ 
          geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
  grid.arrange(g1,g2,g3,g4,ncol=2,nrow=2)
}












