####Rdata to word####
####R to word table
library(stargazer);
data(iris);
head(iris);
stargazer(iris, type = "html", out="D:/star_descriptive.doc");

library(stargazer);
tt <- as.matrix(test.result);
stargazer(tt, title="Table1:Descriptive Statistics",style="asq",notes="Harbin Medical University",type = "html", out="D:/star_descriptive.doc");

mod1 <- lm(Petal.Length ~ Species,data=iris);
mod2 <- lm(Petal.Length ~ Species + Petal.Width,data=iris);
mod3 <- lm(Petal.Length ~ Species + Petal.Width + Species:Petal.Width, data=iris);
stargazer(mod1, mod2, mod3,type="html",out="D:/star_linear.doc");
stargazer(mod1, mod2, mod3,type="html",intercept.bottom = F,intercept.top = T,out="D:/star_linear1.doc");
stargazer(mod1, mod2, mod3,type="html",ci=T,intercept.bottom = F,intercept.top = T,out="D:/star_linear2.doc");
stargazer(mod1, mod2, mod3,type="html",model.names=T,ci=T,intercept.bottom = F,intercept.top = T,out="D:/star_linear3.doc");


##########sjPlot#########
library(sjPlot);
tab_df(iris,title="Table1: Cox jghf dgr dfgr sdgr",show.footnote=T,footnote="Harbin Medical University",CSS=list(css.thead= "border-top:2px solid;",css.lasttablerow = "border-bottom: 2px solid;",css.footnote = "border-top:2px solid;"),file="D:/shijian.doc");
tabl <- tab_df(iris,title="Table1: Cox jghf dgr dfgr sdgr",show.footnote=T,footnote="Harbin Medical University",CSS=list(css.thead= "border-top:2px solid;",css.lasttablerow = "border-bottom: 2px solid;",css.footnote = "border-top:2px solid;"),file="D:/shijian.doc");
cat(tabl$page.content,file="shijianTable.txt");
cat(tabl$page.style,file="shijianTable1CSS.txt");

tabl <- tab_df(test.result,title="Table1: Cox jghf dgr dfgr sdgr",alternate.rows=T,show.footnote=T,footnote="Harbin Medical University",CSS=list(css.thead= "border-top:2px solid;text-align:left;",css.lasttablerow = "border-bottom: 2px solid;",css.footnote = "border-top:2px solid;",css.centeralign = 'text-align: left;',css.firsttablecol = 'font-weight: bold;'),file="D:/shijian.doc");


##########officer#############
library(flextable)
library(officer)
####bg(), fontsize(), italic(), bold(), color(), padding()加样式###
myft <- flextable(head(mtcars)) #输入数据框返回flextable对象
myft <- bold(myft, ~ drat > 3.5, ~ drat, bold = TRUE); #flextable对象中drat列>3.5的加粗
myft <- color(myft, ~ drat > 3.5, ~ drat, color = "red"); #flextable对象中drat列>3.5的变红
doc <- read_docx(); #生成rdoc对象doc
doc <- body_add_flextable(doc, value = myft); #将myft写入rdoc对象doc中
print(doc, target = "test.docx"); 
https://davidgohel.github.io/flextable/articles/overview.html