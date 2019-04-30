library(RCurl)
# 判断url是否存在
url.exists(url="www.baidu.com") # 判断url是否存在
# [1] TRUE
d <- debugGatherer() #收集调试信息
# verbose = TRUE 这时候，d$value()值是会叠加的
tmp <- getURL(url="www.baidu.com", debugfunction = d$update, verbose = TRUE)  
names(d$value())
# [1] "text"       "headerIn"   "headerOut"  "dataIn"     "dataOut"    "sslDataIn"  "sslDataOut"
cat(d$value()[1]) #服务器地址及端口号
cat(d$value()[2]) #服务器返回的头信息
cat(d$value()[3]) #提交给服务器的头信息
d$reset() # 清除d$value()
d$value() # 清除之后全部为空
# text   headerIn  headerOut  dataIn    dataOut  sslDataIn sslDataOut 
# ""         ""         ""         ""         ""         ""         "" 
#2)查看服务器返回的头信息，列表形式
h <- basicHeaderGatherer()
txt <- getURL(url="http://www.baidu.com", headerfunction=h$update)
names(h$value())
cat(h$value())


# getBinaryURL() 下载一个文件
url <- "http://rfunction.com/code/1201/120103.R"
tmp <- getBinaryURL(url) #下载文件
note <- file("D:/120103.R", open = "wb") 
writeBin(tmp, note) #将下载的文件放到指定位置
close(note)



# XML简介
# 缺点：在windows下对中文支持不理想（我在ubuntu下也不理想）
library(XML)
url <- "http://data.earthquake.cn/datashare/datashare_more_quickdata_new.jsp" # 中文界面，抓出来是乱码
url <- "http://219.143.71.11/wdc4seis@bj/earthquakes/csn_quakes_p001.jsp" # 英文界面，抓出来是对的
url <- "http://enhancer.binf.ku.dk/presets/"
wp <- getURL(url)
doc <-htmlParse(wp, asText = TRUE) # 这里切记encoding  
tables <- readHTMLTable(doc, header=F, which = 2) # 选取第二个表
head(tables)



library(RCurl)
library(XML)
webpage <- getURL("http://202.97.205.69/DiseaseEnhancer/")
webpage <- htmlTreeParse(webpage,encoding="GB2312", error=function(...){}, useInternalNodes = TRUE,trim=TRUE)
click <- getNodeSet(webpage,"//a[@href='Quick?query=CARLo-5']")


library(RCurl)
for(i in 1:60){
    Sys.sleep(7200)
	webpage <- getURL("http://202.97.205.69/DiseaseEnhancer/Quick?query=CARLo-5")
}

#R下载命令行
URL <- "http://rfunction.com/code/1202/120222.R"
destfile <- "D:/120222.R"
download.file(URL, destfile)



