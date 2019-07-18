#-----------------------------------
#' 判断芯片表达谱是否需要log2标化
#-----------------------------------
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { 
    ex[which(ex <= 0)] <- NaN
    exprSet <- log2(ex)
    print("log2 transform finished")
}else{
    print("log2 transform not needed")
}
