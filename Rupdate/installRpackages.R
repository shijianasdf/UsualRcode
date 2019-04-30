installRpackages <- function(path){
  #@param path:the path of the packages "C:/Users/lenovo/AppData/Local/Temp/RtmpwlzRsh/downloaded_packages"
  files <- list.files(path,full.names=T);
  for(i in files)
	install.packages(i, repos = NULL, type = "win.binary");
}
installRpackages("C:/Users/lenovo/AppData/Local/Temp/RtmpOiWIiv/downloaded_packages")
installRpackages("C:/Users/lenovo/AppData/Local/Temp/RtmpOiWIiv/downloaded_packages")
installRpackages("C:/Users/lenovo/AppData/Local/Temp/RtmpWCYUg0/downloaded_packages")

#bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GOSemSim", version = "3.8")