library(rvest)
for(i in 1:60){
	Sys.sleep(25200)
	web1 <- read_html("http://biocc.hrbmu.edu.cn/MutExGenome/jsp/Memain!QuickSearch?param=TP53")
	web2 <- read_html("http://bio-bigdata.hrbmu.edu.cn/MutExGenome/jsp/Memain!QuickSearch?param=TP53")
}

