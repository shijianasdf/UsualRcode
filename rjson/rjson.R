library("rjson")
x <- list( alpha = 1:5, beta = "Bravo", 
           gamma = list(a=1:3, b=NULL), 
           delta = c(TRUE, FALSE) )
json <- toJSON( x ) #x是向量或者列表 
fromJSON( json ) #"{\"alpha\":[1,2,3,4,5],\"beta\":\"Bravo\",\"gamma\":{\"a\":[1,2,3],\"b\":null},\"delta\":[true,false]}"
writeLines(json, "D:/R/topic/Meddc/1.MEddc数据/MEddc/网页展示数据/test.json")

y <- 1:10
json <- toJSON(y)
fromJSON( json )
writeLines(json, "D:/R/topic/Meddc/1.MEddc数据/MEddc/网页展示数据/test.json")

x <- list( alpha = 1:5, beta = "Bravo", 
            gamma = list(a=list(1:3), b=NULL), 
            delta = c(TRUE, FALSE) )
json <- toJSON(x) #"{\"alpha\":[1,2,3,4,5],\"beta\":\"Bravo\",\"gamma\":{\"a\":[[1,2,3]],\"b\":null},\"delta\":[true,false]}"



list.list <- list()
for(i in 1:nrow(scat_plot)){
   tempList <- list(name=rownames(scat_plot)[i],data=list(matrix(scat_plot[i,])));
   list.list[[i]] <- tempList
}
json <- toJSON(list.list)
