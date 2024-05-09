Call_scDECODER <- function(dataSet){
  if(length(dataSet$metadata$log.transformed) == 1){
    if(!(dataSet$metadata$log.transformed)){
      tmpEx <- log2(dataSet$ex+1)
    } else{
      tmpEx <- dataSet$ex
    }
  }else{
    tmpEx <- log2(dataSet$ex+1)
  }
  rownames(tmpEx) <- make.names(dataSet$featInfo$SYMBOL, unique = T)
  data <- as.matrix(tmpEx)
  
  # load reference
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  refSet <- read.table("../Rdata/scDECODER/scDECODER_geneWeights.txt", sep="\t")
  geneIDRef <- as.character(refSet[2:nrow(refSet),1])
  geneSigRef <- refSet[2:nrow(refSet),2:ncol(refSet)]
  geneSigRef <- data.frame(sapply(geneSigRef, function(x) as.numeric(as.character(x))))
  names(geneSigRef) <- as.character(unlist(refSet[1,2:ncol(refSet)]))
  
  dataRef <- read.table("../Rdata/scDECODER/sc_exprMatUniq.tsv", sep="\t")
  dataRef <- as.matrix(dataRef[2:nrow(dataRef),2:ncol(dataRef)])
    
  # overlap genes
  geneID <- rownames(data)
  geneOvlp <- geneIDRef[geneIDRef %in% geneID]
  indRef <- match(geneOvlp,geneIDRef)
  dataRef <- dataRef[indRef,]
  geneSigRef <- geneSigRef[indRef,]
  indCur <- match(geneOvlp,geneID)
  data <- data[indCur,]

  # normalize data
  maxRef <- as.numeric(max(dataRef))
  data <- (maxRef/max(data))*data

  # package a nnls function
  callNNLS <- function (A,b){
    res <- nnls::nnls(A,b)$x
    return(res)
  }

  # calculate sample weights
  sampleWeights <- t(apply(data, 2, function(x) callNNLS(as.matrix(geneSigRef),x)))
  colnames(sampleWeights) <- colnames(geneSigRef)

  # normalize
  sampleWeightNorm <- data.frame(sweep(sampleWeights, 1, rowSums(sampleWeights), FUN="/"))
  
  # add parameters
  sampleWeightNorm$bcRatio <- (sampleWeightNorm$Basal_like+0.01)/(sampleWeightNorm$Classical+0.01)
  sampleWeightNorm$bcRatio[sampleWeightNorm$bcRatio > 5] <- 5
  sampleWeightNorm$bcDiff <- sampleWeightNorm$Basal_like-sampleWeightNorm$Classical
  sampleWeightNorm$bcSum <- sampleWeightNorm$Basal_like+sampleWeightNorm$Classical
  
  sampleWeightNorm$miRatio <- (sampleWeightNorm$myCAF+0.01)/(sampleWeightNorm$iCAF+0.01)
  sampleWeightNorm$miRatio[sampleWeightNorm$miRatio > 5] <- 5
  sampleWeightNorm$miDiff <- sampleWeightNorm$myCAF-sampleWeightNorm$iCAF
  sampleWeightNorm$miSum <- sampleWeightNorm$myCAF+sampleWeightNorm$iCAF
  
  dataSet$scDecoderWeights <- sampleWeightNorm
  
  return(dataSet)
}


