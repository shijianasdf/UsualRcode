##### Figure plotting #####
## Figure 1C-D ## 
# load R package
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(patchwork)

# set col
Epithelial.col <- "#FF8214"
Stromal.col <- "#AA40FC"
Myeloid.col <- "#DD4D4D"
Tcells.col <- "#925F54"
Bcells.col <- "#1E77B4"
Mast.col <- "#279E68"
cell.col <- c(Epithelial.col,Stromal.col,Myeloid.col,Tcells.col,Bcells.col,Mast.col)

# load data
dat <- read.delim("meta_annotatedVersion3.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
dat$new.patient <- as.character(dat$patient)
table(dat$new.patient)
dat[which(dat$new.patient == "SSL3"),"new.patient"] <- "HP2"
dat[which(dat$new.patient == "SSL6"),"new.patient"] <- "SSL3"
table(dat$new.patient)
dat$new.patient <- factor(dat$new.patient,levels = rev(c("NC1","NC2","NC3","HP1","HP2","HP3","HP4","SSL1","SSL2","SSL3","SSL4","SSL5","SSLD","TSA1","TSA2","TSA3","TSA4","TSA5")))
dat$annotated <- factor(dat$annotated,levels = c("Epithelial","Stromal cells","Myeloids","T cells","B cells","Mast cells"))

dat.bar <- as.data.frame(table(dat$new.patient,dat$annotated))
dat.bar <- dat.bar%>% group_by(Var1) %>% 
  mutate(Pct = Freq / sum(Freq) * 100) %>% as.data.frame()
colnames(dat.bar) <- c("patient","annotated","count","pct")
dat.bar <- dat.bar[order(dat.bar$patient),]

# left barplot
p.left <- ggplot(dat.bar, aes(x = patient, y = count,fill = annotated)) +
  scale_fill_manual(values  = cell.col) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(name = "",position = "bottom") +
  ylab("Total cell number") + 
  theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.y = element_line(size = 0.6),
        axis.title = element_text(size = 12,colour = "black"),
        axis.text.y = element_text(colour = "black",size = 10),
        axis.title.x = element_text(vjust = -0.3,size = 12),
        axis.text.x = element_text(size = 10, color = "black",angle = 45,hjust = 1),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.title = element_blank())

p.left

p.right <- ggplot(dat.bar, aes(x = patient, y = pct, fill=annotated)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values  = cell.col) +
  scale_x_discrete(name = "") +
  theme_bw() +
  ylab("Total cell proportion (%)") + 
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.y = element_line(size = 0.6),
        axis.title = element_text(size = 12,colour = "black"),
        axis.title.x = element_text(vjust = -0.3,size = 12),
        axis.text.x = element_text(size = 10, color = "black",angle = 45,hjust = 1),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.title = element_blank())
p.right

pal <- p.left + p.right +
  plot_layout(widths = c(6,6), guides = 'collect') & theme(legend.position = 'right',legend.key.size = unit(0.4, 'cm'))
pal
ggsave(filename = "distribution of total cells.pdf", width = 8,height = 5)


# left barplot
dat.bar2 <- as.data.frame(table(dat$group,dat$annotated))
dat.bar2 <- dat.bar2%>% group_by(Var1) %>% 
  mutate(Pct = Freq / sum(Freq) * 100) %>% as.data.frame()
colnames(dat.bar2) <- c("group","annotated","count","pct")
dat.bar2 <- dat.bar2[order(dat.bar2$group),]
dat.bar2$group <- factor(dat.bar2$group,levels = rev(c("TSA","SSLD","SSL","HP","NC")))
p.left <- ggplot(dat.bar2, aes(x = group, y = count,fill = annotated)) +
  scale_fill_manual(values  = cell.col) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(name = "",position = "bottom") +
  ylab("Total cell number") + 
  theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.y = element_line(size = 0.6),
        axis.title = element_text(size = 12,colour = "black"),
        axis.text.y = element_text(colour = "black",size = 10),
        axis.title.x = element_text(vjust = -0.3,size = 12),
        axis.text.x = element_text(size = 10, color = "black",angle = 45,hjust = 1),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.title = element_blank())

p.left

p.right <- ggplot(dat.bar2, aes(x = group, y = pct, fill=annotated)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values  = cell.col) +
  scale_x_discrete(name = "") +
  theme_bw() +
  ylab("Total cell proportion (%)") + 
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.y = element_line(size = 0.6),
        axis.title = element_text(size = 12,colour = "black"),
        axis.title.x = element_text(vjust = -0.3,size = 12),
        axis.text.x = element_text(size = 10, color = "black",angle = 45,hjust = 1),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.title = element_blank())
p.right

pal <- p.left + p.right +
  plot_layout(widths = c(6,6), guides = 'collect') & theme(legend.position = 'right',legend.key.size = unit(0.4, 'cm'))
pal
ggsave(filename = "distribution of total cells2.pdf", width = 5,height = 3)


# use heatmap
dat.bar.wide <- reshape(dat.bar[,c(1,2,4)], idvar = "patient", timevar = "annotated", direction = "wide")
rownames(dat.bar.wide) <- dat.bar.wide$patient; dat.bar.wide <- dat.bar.wide[,-1]
colnames(dat.bar.wide) <- gsub("pct.","",colnames(dat.bar.wide),fixed = T)
dat.bar.wide <- dat.bar.wide/100
rowSums(dat.bar.wide)

hm.dat <- data.frame(row.names = c("TSA1","TSA2","TSA3","TSA4","TSA5",
                                   "SSLD","SSL1","SSL2","SSL3","SSL4","SSL5",
                                   "HP1","HP2","HP3","HP4",
                                   "NC1","NC2","NC3"))
hm.dat <- matrix(0,nrow = 18,ncol = 6, byrow = T,
                 dimnames = list(c("TSA1","TSA2","TSA3","TSA4","TSA5",
                                   "SSLD","SSL1","SSL2","SSL3","SSL4","SSL5",
                                   "HP1","HP2","HP3","HP4",
                                   "NC1","NC2","NC3"),
                                 c("Rectum","Sigmoid","Transverse colon","Hepatic flexure","Ascending colon","Cecum")))
hm.dat <- as.data.frame(hm.dat)
hm.dat["TSA1","Rectum"] <- 1
hm.dat[c("TSA2","TSA3"),"Sigmoid"] <- 1
hm.dat[c("TSA4","TSA5","SSLD","SSL4","HP3","NC1","NC2","NC3"),"Transverse colon"] <- 1
hm.dat[c("SSL1","SSL3","HP4"),"Ascending colon"] <- 1
hm.dat[c("SSL2","HP1","HP2"),"Cecum"] <- 1
hm.dat[c("SSL5"),"Hepatic flexure"] <- 1
dat.bar.wide <- dat.bar.wide[rownames(hm.dat),]

count <- reshape(dat.bar[,c(1,2,3)], idvar = "patient", timevar = "annotated", direction = "wide")
rownames(count) <- count$patient; count <- count[,-1]
colnames(count) <- gsub("count.","",colnames(count),fixed = T)
rowSums(count)

annRow <- data.frame(row.names = rownames(hm.dat),
                     group = factor(rep(c("TSA","SSLD","SSL","HP","NC"),c(5,1,5,4,3)),levels = c("TSA","SSLD","SSL","HP","NC")),
                     stringsAsFactors = F)
annRow$color <- annColors$group[annRow$group]
annRow$count <- rowSums(count[rownames(annRow),])
annColors <- list(group = c("TSA" = "#D6604D",
                            "SSLD" = "#689E45",
                            "SSL" = "#5A8CA8",
                            "HP" = "#93C1B6",
                            "NC" = "#DFC27D"))
group.col <- c("#D6604D","#689E45","#5A8CA8","#93C1B6","#DFC27D"); names(group.col) <- c("TSA","SSLD","SSL","HP","NC")
sample.col <- annRow$color; names(sample.col) <- annRow$group

hm <- Heatmap(as.matrix(hm.dat),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_names_side = "top",
        show_heatmap_legend = FALSE,
        border = FALSE,
        left_annotation = HeatmapAnnotation(which = "row",
                                            Group = anno_simple(factor(annRow$group,levels = c("TSA","SSLD","SSL","HP","NC")),
                                                                col = sample.col),
                                            Proportion = anno_barplot(as.matrix(dat.bar.wide),
                                                                      which = "row",
                                                                      border = FALSE,
                                                                      gp = gpar(fill = cell.col,border=NA,lty="blank"), 
                                                                      bar_width = 0.85,
                                                                      axis_param = list(side = "top",at = c(0,0.25,0.5,0.75,1),labels = c(0,25,50,75,100)),
                                                                      width = unit(4, "cm"),
                                                                      height = unit(1, "cm")),
                                            show_annotation_name = FALSE),
        # left_annotation = rowAnnotation(Group = factor(annRow$group,levels = c("TSA","SSLD","SSL","HP","NC")),
        #                                 col = list(Group = sample.col),
        #                                 show_annotation_name = FALSE,
        #                                 annotation_label = "Patient category"),
        right_annotation = HeatmapAnnotation(which = "row",
                                             Count = anno_points(annRow$count,
                                                                 pch = 19, 
                                                                 col = "black",
                                                                 border = FALSE,
                                                                 width = unit(1,"cm"),
                                                                 axis_param = list(side = "top",at = c(2000,4000,6000),labels = c(2000,4000,6000))),
                                             show_annotation_name = FALSE),
        
        rect_gp = gpar(col = "white"),
        width = ncol(hm.dat)*unit(5, "mm"), 
        height = nrow(hm.dat)*unit(5, "mm"),
        row_split = annRow$group,
        gap = unit(2,"mm"),
        col = c("grey90","grey30"))

lgd_group = Legend(labels = c("TSA","SSLD","SSL","HP","NC"), title = "Patient category", legend_gp = gpar(fill = unique(sample.col)))
lgd_proportion = Legend(labels = colnames(dat.bar.wide), title = "Cell type", legend_gp = gpar(fill = cell.col))
pdf(file = "landscape of cell proportion in samples.pdf",width = 6,height = 6)
draw(hm, annotation_legend_list = list(lgd_group, lgd_proportion))
invisible(dev.off())



## Figure 1E ## 
load("seurat_full.Rdat")
plot3<-FeaturePlot(seurat, features = c("EPCAM", "CD3D","CD14",'CD163',"CD79A", "MS4A1", "COL1A1",'PECAM1'),ncol =4,
                   cols=c('#D3D3D3','#FF8D28'))#Landscape
ggsave("serrated adenoma/Featureplot 20211213-1.pdf", plot = plot3, width = 15, height = 6)

plot4<-FeaturePlot(seurat, features = c("EPCAM", "CD3D","CD14",'CD163',"CD79A", "MS4A1", "COL1A1",'PECAM1'),ncol =4,
                   cols=c('#D3D3D3','#8C564B'))#Landscape
ggsave("serrated adenoma/Featureplot 20211213-2.pdf", plot = plot4, width = 15, height = 6)

plot5<-FeaturePlot(seurat, features = c("EPCAM", "CD3D","CD14",'CD163',"CD79A", "MS4A1", "COL1A1",'PECAM1'),ncol =4,
                   cols=c('#D3D3D3','#DD4D4D'))#Landscape
ggsave("serrated adenoma/Featureplot 20211213-3.pdf", plot = plot5, width = 15, height = 6)

plot6<-FeaturePlot(seurat, features = c("EPCAM", "CD3D","CD14",'CD163',"CD79A", "MS4A1", "COL1A1",'PECAM1'),ncol =4,
                   cols=c('#D3D3D3','#2178B5'))#Landscape
ggsave("serrated adenoma/Featureplot 20211213-4.pdf", plot = plot6, width = 15, height = 6)

plot7<-FeaturePlot(seurat, features = c("EPCAM", "CD3D","CD14",'CD163',"CD79A", "MS4A1", "COL1A1",'PECAM1'),ncol =4,
                   cols=c('#D3D3D3','#AC43FC'))#Landscape
ggsave("serrated adenoma/Featureplot 20211213-5.pdf", plot = plot7, width = 15, height = 6)

