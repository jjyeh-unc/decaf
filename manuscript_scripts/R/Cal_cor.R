Cal_cor <- function(vector, centroidMat, colStart, colEnd) {
  idxMax <- colStart
  rMax <- -999
  idxMax <- colStart
  for (i in colStart:colEnd) {
    rTmp <- cor(vector, centroidMat[,i], method = "pearson")
    if(rTmp>rMax) { 
      rMax <- rTmp ;
      idxMax <- i;}
  }
  res <- data.frame(idx = idxMax,
                    r = rMax,
                    stringsAsFactors = FALSE)
  return(res)
}