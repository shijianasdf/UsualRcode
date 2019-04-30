AutoCreateDir <- function(rootPath,Dirfile){
	paths <- file.path(rootPath,Dirfile$firstDir,Dirfile$secondDir)
	for(i in 1:length(paths)){
		if(!file.exists(paths[i])){
			dir.create(paths[i],recursive=T)
		}
	}
}
