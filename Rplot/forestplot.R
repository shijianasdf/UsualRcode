#############################################
# Simple examples of how to do a forestplot #
#############################################
library(forestplot)
ask <- par(ask=TRUE)

# A basic example, create some fake data
row_names <- list(list("test = 1", expression(test >= 2)))
test_data <- data.frame(coef=c(1.59, 1.24),
                        low=c(1.4, 0.78),
                        high=c(1.8, 1.55))
forestplot(row_names,
           test_data$coef,
           test_data$low,
           test_data$high,
           zero = 1,
           cex  = 2,
           lineheight = "auto",
           xlab = "Lab axis txt")

# Print two plots side by side using the grid
# package's layout option for viewports
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.col = 1))
forestplot(row_names,
           test_data$coef,
           test_data$low,
           test_data$high,
           zero = 1,
           cex  = 2,
           lineheight = "auto",
           xlab = "Lab axis txt",
           new_page = FALSE)
popViewport()
pushViewport(viewport(layout.pos.col = 2))
forestplot(row_names,
           test_data$coef,
           test_data$low,
           test_data$high,
           zero = 1,
           cex  = 2,
           lineheight = "auto",
           xlab = "Lab axis txt",
           new_page = FALSE)
popViewport(2)


# An advanced test
test_data <- data.frame(coef1=c(1, 1.59, 1.3, 1.24),
                        coef2=c(1, 1.7, 1.4, 1.04),
                        low1=c(1, 1.3, 1.1, 0.99),
                        low2=c(1, 1.6, 1.2, 0.7),
                        high1=c(1, 1.94, 1.6, 1.55),
                        high2=c(1, 1.8, 1.55, 1.33))

col_no <- grep("coef", colnames(test_data))
row_names <- list(
  list("Category 1", "Category 2", "Category 3", expression(Category >= 4)),
  list("ref",
       substitute(expression(bar(x) == val),
                  list(val = round(rowMeans(test_data[2, col_no]), 2))),
       substitute(expression(bar(x) == val),
                  list(val = round(rowMeans(test_data[3, col_no]), 2))),
       substitute(expression(bar(x) == val),
                  list(val = round(rowMeans(test_data[4, col_no]), 2))))
)

coef <- with(test_data, cbind(coef1, coef2))
low <- with(test_data, cbind(low1, low2))
high <- with(test_data, cbind(high1, high2))
forestplot(row_names, coef, low, high,
           title="Cool study",
           zero = c(0.98, 1.02),
           grid = structure(c(2^-.5, 2^.5), gp = gpar(col = "steelblue", lty=2)),
           boxsize=0.25,
           col=fpColors(box=c("royalblue", "gold"),
                        line=c("darkblue", "orange"),
                        summary=c("darkblue", "red")),
           xlab="The estimates",
           new_page = TRUE,
           legend=c("Treatment", "Placebo"),
           legend_args = fpLegend(pos = list("topright"),
                                  title="Group",
                                  r = unit(.1, "snpc"),
                                  gp = gpar(col="#CCCCCC", lwd=1.5)))

# An example of how the exponential works
test_data <- data.frame(coef=c(2.45, 0.43),
                        low=c(1.5, 0.25),
                        high=c(4, 0.75),
                        boxsize=c(0.5, 0.5))
row_names <- cbind(c("Name", "Variable A", "Variable B"),
                   c("HR", test_data$coef))
test_data <- rbind(rep(NA, 3), test_data)

forestplot(labeltext = row_names,
           test_data[,c("coef", "low", "high")],
           is.summary=c(TRUE, FALSE, FALSE),
           boxsize   = test_data$boxsize,
           zero      = 1,
           xlog      = TRUE,
           col = fpColors(lines="red", box="darkred"))


par(ask=ask)
# See vignette for a more detailed description
# vignette("forestplot",  package="forestplot")


workdir <- "C:\\Path\\To\\Relevant\\Directory"
datafile <- file.path(workdir,"ForestPlotData.csv")
data <- read.csv(datafile, stringsAsFactors=FALSE)

## Labels defining subgroups are a little indented!
subgps <- c(4,5,8,9,12,13,16,17,20,21,24,25,28,29,32,33)
data$Variable[subgps] <- paste("  ",data$Variable[subgps]) 

## Combine the count and percent column
np <- ifelse(!is.na(data$Count), paste(data$Count," (",data$Percent,")",sep=""), NA)

## The rest of the columns in the table. 
tabletext <- cbind(c("Subgroup","\n",data$Variable), 
                   c("No. of Patients (%)","\n",np), 
                   c("4-Yr Cum. Event Rate\n PCI","\n",data$PCI.Group), 
                   c("4-Yr Cum. Event Rate\n Medical Therapy","\n",data$Medical.Therapy.Group), 
                   c("P Value","\n",data$P.Value))

library(forestplot)
png(file.path(workdir,"Figures\\Forestplot.png"),width=960, height=640)
forestplot(labeltext=tabletext, graph.pos=3, 
           mean=c(NA,NA,data$Point.Estimate), 
           lower=c(NA,NA,data$Low), upper=c(NA,NA,data$High),
           title="Hazard Ratio",
           xlab="     <---PCI Better---    ---Medical Therapy Better--->",
           hrzl_lines=list("3" = gpar(lwd=1, col="#99999922"), 
                           "7" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
                           "15" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
                           "23" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
                           "31" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922")),
           txt_gp=fpTxtGp(label=gpar(cex=1.25),
                          ticks=gpar(cex=1.1),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           col=fpColors(box="black", lines="black", zero = "gray50"),
           zero=1, cex=0.9, lineheight = "auto", boxsize=0.5, colgap=unit(6,"mm"),
           lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.4)
dev.off()
