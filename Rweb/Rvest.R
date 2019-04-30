library(rvest)
rootpath <- "http://tfbsbank.co.uk" 
TFBSbank <- "D:/R/topic/Network/数据/TFBSbank"#下载数据的根目录
web <- read_html("http://tfbsbank.co.uk/view.php?species=Homo%20sapiens")#首页面
webtable <- html_table(web)[[1]]
path <- paste(webtable$TF,webtable$ID,webtable$Stage_tissue,sep="_")
dirpath <- file.path(TFBSbank,path) #下载数据的文件夹路径
for(i in 1:length(dirpath)){
	if(!file.exists(dirpath[i]))
	dir.create(dirpath[i],showWarnings=TRUE,recursive=TRUE,mode="775") 
}
text <- web %>% html_nodes("table#customers") %>% html_nodes("tr>td>a[href]") %>% html_attr("href") #创建文件夹
jumpToDownloadPage <- file.path(rootpath,text) #
for(i in 1:length(jumpToDownloadPage)){
	web <- read_html(jumpToDownloadPage[i])
	downloadUrl <- web %>% html_nodes("div#introduction") %>% html_nodes("ul>li>a[href]") %>% html_attr("href")
	destfile <- file.path(dirpath[i],basename(jumpToDownloadPage[i]))
    download.file(downloadUrl,destfile)
	#tmp <- getBinaryURL(downloadUrl,ssl.verifypeer=FALSE) #下载文件
    #note <- file(file.path(dirpath[i],basename(jumpToDownloadPage[i])), open = "wb") 
    #writeBin(tmp, note) #将下载的文件放到指定位置
    #close(note)
}

library(rvest)
web <- read_html("http://firebrowse.org/")
web %>% html_nodes("select#cohort-name-select") %>% html_nodes("option")



