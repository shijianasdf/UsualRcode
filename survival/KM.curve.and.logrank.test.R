#函数一 grouping.bymutation.singleMarker
#' @description 依据单个基因是否突变对样本进行分类
#' @param markers marker基因的名字
#' @param clinical.data  干净的临床数据，每行代表一个样本，每列代表样本的不同临床信息
#' @param mutation.data  干净的突变数据，每行代表一个突变，每列代表突变的相关信息
#' @return 返回的是带有标签的临床数据，最后一列sample.label为分类标签，是通过这个函数产生的，必须为factor
grouping.bymutation.singleMarker <- function(markers, 
                                  clinical.data, 
                                  mutation.data
)
{
  #STEP 1:处理markers，要求markers必须有效
  markers <- intersect(markers ,unique(mutation.data$Gene_Symbol));
  if(length(markers) == 0)
  {
    stop("所选则的marker并未出现在突变谱中！")
  }
  
  #STEP 2:处理样本，保留在突变数据和临床数据中共同出现的样本
  common.sample <- intersect(unique(mutation.data$Patient_ID), clinical.data$Patient_ID );
  clinical.data <- clinical.data[clinical.data$Patient_ID %in% common.sample, ];
  mutation.data <- mutation.data[mutation.data$Patient_ID %in% common.sample, ];
  
  #STEP 3:对样本进行分类，根据是否突变将样本分为突变型和野生型
  label <- c(paste(markers[1], "WT",sep = " "),paste(markers[1], "mut",sep = " "));
  mut.sample <- unique(mutation.data[mutation.data$Gene_Symbol %in% markers,"Patient_ID"]);
  clinical.data$'sample.label' <- factor(label[as.numeric(clinical.data$Patient_ID %in% mut.sample)+1], levels = label);
  return(clinical.data);
}


#函数二 grouping.bymutation.twoMarkers
#' @description 依据两个marker基因是否突变对样本进行分类
#' @param markers marker基因的名字
#' @param clinical.data  干净的临床数据，每行代表一个样本，每列代表样本的不同临床信息，至少包含三列：Patient_ID表示临床样本ID；event表示结局事件；time表示样本的结局事件时间
#' @param mutation.data  干净的突变数据，每行代表一个突变，每列代表突变的相关信息，至少包含两列：Gene_Symbol表示发生突变的基因名；Patient_ID表示对应突变所在的样本ID（必须与临床信息的样本ID一致）
#' @return 返回的是带有标签的临床数据，最后一列sample.label为分类标签，是通过这个函数产生的，必须为factor
grouping.bymutation.twoMarkers <- function(markers, 
                                           mutation.data, 
                                           clinical.data
)
{
  
    #STEP 1:处理markers，要求markers必须有效
    markers <- intersect(markers ,unique(mutation.data$Gene_Symbol));
    if(length(markers) == 0)
    {
        stop("所选则的marker并未出现在突变谱中！")
    }
  
  
    #STEP 2:处理样本，保留在突变数据和临床数据中共同出现的样本
    common.sample <- intersect(unique(mutation.data$Patient_ID), clinical.data$Patient_ID );
    clinical.data <- clinical.data[clinical.data$Patient_ID %in% common.sample, ];
    mutation.data <- mutation.data[mutation.data$Patient_ID %in% common.sample, ];
  
  
    #STEP 3: 依据特定的方法将样本分类，分配分类标签
  
    mut.sample1 <- unique(mutation.data[mutation.data$Gene_Symbol %in% markers[1],"Patient_ID"]);
    mut.sample2 <- unique(mutation.data[mutation.data$Gene_Symbol %in% markers[2],"Patient_ID"]);
    
    label1 <- c(paste(markers[1],"WT",sep = ""),paste(markers[1],"mut",sep = ""));
    label2 <- c(paste(markers[2],"WT",sep = ""),paste(markers[2],"mut",sep = ""));
    
    sample.label1 <- factor(label1[as.numeric(clinical.data$Patient_ID %in% mut.sample1)+1], levels = label1);
    sample.label2 <- factor(label2[as.numeric(clinical.data$Patient_ID %in% mut.sample2)+1], levels = label2);
    sample.label <- paste(sample.label1, sample.label2, sep="&");  
  
    #STEP 4：更换样本标签格式为常用的类型，并添加标签到clinical.data中，记作sample.label 

    label.level <- c("WT", paste("only", markers[1],"mut"), paste("only", markers[2],"mut"), paste(markers[1],"and", markers[2],"mut"));
    
    sample.label[which(sample.label == paste(markers[1],"WT","&",markers[2],"WT",sep = ""))] <- label.level[1]
    sample.label[which(sample.label == paste(markers[1],"WT","&",markers[2],"mut",sep = ""))] <- label.level[3];
    sample.label[which(sample.label == paste(markers[1],"mut","&",markers[2],"WT",sep = ""))] <- label.level[2];
    sample.label[which(sample.label == paste(markers[1],"mut","&",markers[2],"mut",sep = ""))] <- label.level[4];
    
    clinical.data$sample.label <- factor(sample.label, levels = label.level);
  
    return(clinical.data)
}
#
#'-------------------------------------------------
#' @author jian shi
#' @description 对每个病人的打上对应基因的CCF值
#' @param  marker 基因名字
#' @param  clinical.data 临床数据
#' @param  clone.mutation.data 克隆突变数据
#' @example lable.byCCF.Markers("APC",)
#'-------------------------------------------------

lable.byCCF.Markers <- function(marker.symbol,clinical.data,clonality.data){
  #STEP 1: 处理样本，保留在突变数据和临床数据中共同出现的样本
  common.sample <- intersect(unique(clonality.data$Patient_ID), clinical.data$Patient_ID);
  clinical.data <- clinical.data[clinical.data$Patient_ID %in% common.sample, ];
  clonality.data <- clonality.data[clonality.data$Patient_ID %in% common.sample, ];
  
  #STEP 2:识别marker基因是否出现在克隆性数据中
  markers <- intersect(marker.symbol,unique(clonality.data$Gene_Symbol));
  if(length(markers) == 0)
  {
    sample.label <- rep(paste(marker.symbol,"WT"),length(common.sample));
    return(cbind.data.frame(clinical.data,sample.label));
  }  
  
  #STEP 3: 得到克隆、亚克隆样本
  clonality.cellmarker.data <- unique(clonality.data[clonality.data$Gene_Symbol %in% marker.symbol,])   
  clone.sample <- unique(clonality.cellmarker.data[which(clonality.cellmarker.data$"CI95.timing" == "Clonal"),1])
  subclone.sample <- unique(clonality.cellmarker.data[which(clonality.cellmarker.data$"CI95.timing" == "Subclonal"),1])
  #####找到对于某个基因即是克隆又是亚克隆的样本,在我们分析中将这样的克隆亚克隆样本去掉,并且给样本打上标签
  temp.sample <- intersect(clone.sample,subclone.sample)
  if(length(temp.sample)==0){
    label1 <- c(paste(markers, "WT",sep = " "),paste(markers, "Clonal",sep = " "));
    label2 <- c(paste(markers, "WT",sep = " "),paste(markers, "Subclonal",sep = " "));
    level <- c(paste(markers, "WT",sep = " "),paste(markers, "Clonal",sep = " "),paste(markers, "Subclonal",sep = " ")); #设置因子顺序，对应画图颜色
    sample.label1 <- label1[as.numeric(clinical.data$Patient_ID %in% clone.sample)+1];
    sample.label2 <- label2[as.numeric(clinical.data$Patient_ID %in% subclone.sample)+1];
    pos2 <- grep("Subclonal",sample.label2);
    sample.label1[pos2] <- paste(markers,"Subclonal",sep=" ");
    clinical.data$'sample.label' <- factor(sample.label1,levels=level); #数据框新建一列	 
  }else{
    label1 <- c(paste(markers, "WT",sep = " "),paste(markers, "Clonal",sep = " "));
    label2 <- c(paste(markers, "WT",sep = " "),paste(markers, "Subclonal",sep = " "));
    level <- c(paste(markers, "WT",sep = " "),paste(markers, "Clonal",sep = " "),paste(markers, "Subclonal",sep = " ")); #设置因子顺序，对应画图颜色	 
    sample.label1 <- label1[as.numeric(clinical.data$Patient_ID %in% clone.sample)+1];
    sample.label2 <- label2[as.numeric(clinical.data$Patient_ID %in% subclone.sample)+1];
    pos1 <- grep("Clonal",sample.label1);
    pos2 <- grep("Subclonal",sample.label2);
    Clonal.Subclonal.pos <- intersect(pos1,pos2);
    sample.label1[pos2] <- paste(markers,"Subclonal",sep=" ");
    sample.label1[Clonal.Subclonal.pos] <- paste(markers,"WT",sep = " ");
    clinical.data$'sample.label' <- factor(sample.label1,levels=level); #数据框新建一列	 
  }              
  return(clinical.data)
}

#函数三 plot.surv
#' @description 绘制KM曲线并且进行logrank检验
#' @param clinical.data 具有分类标签的临床数据信息，至少包含四列：Patient_ID表示临床样本ID；event表示样本事件结局(数值型或逻辑型)；time表示样本事件的时间（数值型）；最后一列是通过上面的分类函数产生的分类标签sample.label，必须是factor
#' @param upper.time: 数值型，生存时间上限，必须与clinical.data中的单位一致，默认为NULL；如果使用，超过年限的样本将被去掉
#' @param xscale: 字符型，允许的选项包括"d_m"，"d_y"，"m_d"，"m_y"，"y_d"和"y_m"，其中d =天，m =月和y =年。 例如，xscale =“d_m”会将x轴单位从几天转换为几个月
#' @param unit.xlabel: x轴label展示的时间单位，四种选择：c("year", "month", "week", "day"),分别代表年，月，周，天，必须与最终图中展示使用的时间单位一致
#' @param surv.median.line: 在中位生存期时绘制水平/垂直线的字符向量，允许的值包括c（"none"，"hv"，"h"，"v"）中的一个，v：垂直，h：水平
#' @param risk.table: 逻辑向量，是否绘制risk.table，默认为TRUE
#' @param pval: 逻辑向量，是否给出log rank p值，默认为TRUE
#' @param survival.event: 结局事件的类型，允许的值包括c(“Overall Survival”, “Progress Free Survival”)等
#' @param main: 主标题的名字
#' @param inputFilePath: KMplot存储路径
#' @param picture.name: 存储的KMplot名字
#' @return 中位生存期，1,3,5年生存率，KM曲线图片

plot.surv <- function(clinical.data,
                      upper.time=NULL,
                      xscale = c( "d_m", "d_y", "m_d", "m_y", "y_d", "y_m"),#用于转换x轴坐标刻度
                      unit.xlabel = c("year", "month", "week", "day"),#所用生存时间的单位,
                      surv.median.line = c("none", "hv", "h", "v"),#是否画出中位生存时间，默认不给出
                      risk.table = c(TRUE, FALSE),#是否显示risk table
                      pval = c(TRUE, FALSE),#是否给出log-rank的p值
                      conf.int = c(FALSE, TRUE),#是否画出置信区间
                      main=NULL,#主标题名字
                      survival.event=c("Overall Survival","Progress Free Survival"),#事件类型
					            inputFilePath = NULL, #KMplot存储路径
					            picture.name = NULL  #存储的KMplot名字
)
{
    options(stringsAsFactors = FALSE);
    #载入生存分析包
    require(survival);
	  # ggplot2 survival plot ggsurvplot and pairwise_survdiff
    require(survminer);
    require(RColorBrewer);
	  # Provides a number of user-level functions to work with "grid" graphics, notably to arrange multiple grid-based plots on a page, and draw tables.
    require(gridExtra);
    
    survival.event <- survival.event[1];
    unit.xlabel <- unit.xlabel[1];
    
    #去除生存时间超过upper.time的样本
    if(!is.null(upper.time))
    {
        clinical.data <- clinical.data[clinical.data$time <= upper.time,]
    }
    
    #选择日期格式 
    xSL <- data.frame(xScale=c(1,7,30,365.25),xLab=c("Days","Weeks","Months","Years"), stringsAsFactors=FALSE)
    switch(unit.xlabel, year={xScale <- 365.25;}, month={xScale <- 30;}, week={xScale <- 7;}, day={xScale <- 1})
    xLab <- xSL[which(xSL[,1]==xScale),2];
    
    #构造颜色
    t.name <- levels(clinical.data$sample.label);
    colors <- c("#808080","#EA4335","#4285F4","#FBBC05","#34A853","#000000");#顺序：灰，红，蓝，黄，绿，黑
    t.col <- colors[1:length(t.name)]; #palette = c("#E7B800", "#2E9FDF")
    
    #构造生存对象Surv(time, event)并且依据sample.label分层绘制KM生存曲线
    km.curves <- survfit(Surv(time, event)~sample.label, data=clinical.data);
    #km.curves <- survfit(Surv(time, event)~sample.label, data=OS.byDriverGene.ClonalSubclonal[[1]]);
	  
    #生成log-rank检验的p值,保留两位有效数字
	  pval <- round(surv_pvalue(km.curves, clinical.data)$pval,digits=4);
    
	  #计算中位生存期以及CI
    median.os <- surv_median(fit=km.curves);
	
	  #计算5年生存率
    if(unit.xlabel == "month")
	    survival_rate <- summary(km.curves,times=c(12,36,60)); 
    if(unit.xlabel == "year")
      survival_rate <- summary(km.curves,times=c(1,3,5));
    if(unit.xlabel == "week")
      survival_rate <- summary(km.curves,times=c(104.35,156.53,260.89));
    if(unit.xlabel == "day")
      survival_rate <- summary(km.curves,times=c(730.5,1095.75,1826.25));

    #返回中位生存期和5年生存率
    if(pval <= 0.05){
      result <- list(median.os=median.os,survival_rate=survival_rate);
    }else{
      result <- NULL;
    }
      
	  #构造生存图像中图例显示文字
    legend.content <- substr(names(km.curves$strata),start = 14,stop = 1000);
    
  	#设置生存图片地址及名字
  	if(is.null(inputFilePath)){
  	   inputFilePath <- getwd();
  	}
  	if(!file.exists(inputFilePath)){
  	   dir.create(inputFilePath,recursive=T);
  	}
  	if(is.null(picture.name)){
  		picture.name <- "default";
  	}
  	pdf(paste0(inputFilePath,picture.name,"_",pval,".pdf"));
	
    #ggsurvplot绘制生存图像
    ggsurv <- ggsurvplot(                         #注意此时并没有画图，只是存在ggsurv变量中
                         km.curves,               # survfit object with calculated statistics.
                         data = clinical.data,             # data used to fit survival curves.
                         palette = t.col,
                         
                         #图的主题构架
                         risk.table = risk.table[1],       # 是否展示风险table.
                         pval = pval[1],             # 是否展示log-rank检验的p值.
                         surv.median.line = surv.median.line[1],  # 是否展示中位生存期.surv_median
						             conf.int = conf.int[1], #是否画出置信区间
                         title = main,     #主标题名字
                         font.main = 15,       #主标题字体大小              
                         xlab = paste("Time","in",xLab,sep = " "),   # customize X axis label.
                         ylab = survival.event,   # customize Y axis label.
                         
                         #图例设置
                         legend.title = "", #图例标题，一般不用，设为空
                         legend.labs = legend.content, #图例文字描述
                         #legend = c(0.8,0.9), #图例的位置，取值在【0,1】之间
                         font.legend = 9,     #图例字体大小
                         
                         #risk table设置
                         tables.theme = theme_cleantable(),#table主题
                         risk.table.title = "No. at risk:",#table标题
                         risk.table.y.text.col = T, # 使用颜色代替Y轴文字
                         risk.table.y.text = FALSE, # Y轴不使用文字注释
                         tables.height = 0.15,      # table的高度
                         risk.table.fontsize = 3,    # risk table内文字的大小
						             ggtheme = theme_bw() # Change ggplot2 theme
                        );
    ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"));#KM曲线图的标题居中
    ggsurv$table <- ggsurv$table + theme(plot.title = element_text(hjust = -0.04), plot.margin = unit(c(5.5, 5.5, 5.5, 50), "points"));#风险table图的标题居中
    
    #判断分类的类数，如果只有两类，就不必计算两两之间的log rank p值
    if(length(t.name) > 2)
    {
    
    #计算pairwise的log rank的p值
    res <- pairwise_survdiff(Surv(time, event)~sample.label, data=clinical.data);
    pairwise.pvalue <- round(res$p.value, digits = 4);
    pairwise.pvalue[which(pairwise.pvalue < 0.0001)] <- "<0.0001";
    pairwise.pvalue[is.na(pairwise.pvalue)] <- "-";
    
    #添加表格
    tt <- ttheme_minimal(
                         core = list(fg_params = list(col = "black"),bg_params = list(fill = NA, col = "black")),
                         colhead = list(fg_params = list(col = NA),bg_params = list(fill = t.col, col = "black")),
                         rowhead = list(fg_params = list(col = NA, hjust = 1),bg_params = list(fill = c("white",t.col[-1]), col = "black"))
                        );
    pairwise.table <- tableGrob(pairwise.pvalue, theme = tt);
    print(ggarrange(ggarrange(ggsurv$plot, ggsurv$table,nrow=2,heights=c(2,0.5)), pairwise.table, nrow=2,heights = c(2,0.5),labels = c("","p from pairwise comparisons"), hjust = 0, font.label = list(size = 15, face = "plain"))); #最终的绘图函数ggarange(Arrange multiple ggplots on the same page)，包括3部分ggsurv$plot(KM生存曲线)，ggsurv$table(风险table)，成对比较table，print是为了让循环画图可以进行
    }
    else
    {
    print(ggsurv)
    }
	dev.off();
	#返回中位生存期和5年生存率
	return(result);
	#ggsave("shijian.pdf", device="pdf",dpi=300, plot=a, path="D:/");也可以用ggsave保存图片 
}
