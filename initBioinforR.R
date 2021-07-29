#' @title easy report system in function building
#' @description easy report system in function building
#' @param levels an integer >= 1
#' @param type one of "cat" and "message"
#' @param ... one or multiple characters
#' @return a verbose report
#' @seealso \code{\link[base]{cat}};\code{\link[base]{message}}
#' @author shi jian<\email{280462610@qq.com}>
#' @examples
#' LuckyVerbose("AVC")
#' LuckyVerbose("AVC",levels = 1)
#' LuckyVerbose("AVC",levels = 3)
#' LuckyVerbose("AVC",levels = 3,type="message")
#' @export
LuckyVerbose <- function(...,levels = 1,type = NULL){
  ## Verbose type
  if(is.null(type)){
    if(levels == 1){
      type <-  "message"
    } else {
      type <-  "cat"
    }
  }

  ## level symbol
  if(levels > 1){
    s1 <- paste(rep(" ",(levels-1)),collapse = "")
    s2 <- paste(rep("o",(levels-1)),collapse = "")
    ls <- paste(s1,s2,collapse = "")
  } else {
    ls <- ""
  }

  ## do Verbose
  if(type == "message"){
    return(base::message(ls," ",...))
  } else {
    if(type == "cat"){ return(base::cat(ls,...,"\n")) } else {
      print("Input right type.")
    }
  }

}
## QuickStart help run series R script automatically
# script.path# the path of script
# data.path#the path of data
# result.path#the path of result
#' @export
QuickStart <- function(script.path,
                       data.path,
                       result.path){
  ## create result path
  if(!file.exists(result.path)){
    LuckyVerbose(paste0("Create ",result.path))
    dir.create(result.path,recursive = T)
  }
  result.path <<- result.path
  data.path <<- data.path
  ## run scripts
  s <- list.files(path = script.path,pattern = ".R",full.names = T)
  for(i in s) source(i,encoding = "UTF-8")
  
}
QuickStart("/Users/biofly/Rcode/BioinforR/R","/Users/biofly/Rcode/BioinforR/data","/Users/biofly/Rcode/BioinforR/result")

#' @export
## load1 help load local .rda files by a pattern
load1 <- function(pattern,path=".",envir=parent.frame()){
  for(i in pattern){
    file = list.files(path = path,pattern = i,full.names = T)
    if(length(file) > 1){
      print("Attention!Multiple files.")
    }
    for(j in file){
      load(j,envir = envir)
    }
  }
}
load1(".rda",path="/Users/biofly/Rcode/BioinforR/data",envir=parent.frame())
