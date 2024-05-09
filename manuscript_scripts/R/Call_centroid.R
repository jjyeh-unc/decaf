Call_centroid <- function(dataSet, tmpSchema, sampSub, lowVarFlt, distance="euclidean", Rversion = "R-4.X.X") {
  dataSet$tmp <- list()
  dataSet$tmp$mainLab <- tmpSchema
  dataSet$tmp$schema <- tmpSchema
  dataSet$tmp$sampSub <- sampSub
  
  ex <- log2(dataSet$ex+1)
  
  # get genes
  i <- which(schemaList %in% tmpSchema)
  centroidMat <- subtypeGeneList[[i]]
  
  # filter undetected genes
  genesTmp <- dataSet$featInfo$SYMBOL[which(dataSet$featInfo$SYMBOL %in% centroidMat$geneSymbol)]
  centroidMat <- centroidMat[which(centroidMat$geneSymbol %in% genesTmp),]
  
  # match genes in centroid
  idxGene <- match(centroidMat$geneSymbol,dataSet$featInfo$SYMBOL)
  ex <- ex[idxGene,]
  
  # filter samples
  ex <- ex[, sampSub]
  dataSet$tmp$ex <- ex
  
  # calculate corr
  clusterRes <- do.call(rbind, 
                        apply(ex, MARGIN=2, 
                              FUN=function(x) Cal_cor(x, centroidMat, 2, 6)))
  clusterRes$subtypeCall <- names(centroidMat)[clusterRes$idx]
  
  # get orders of sample
  sampSrtAll <- c()
  for (j in 2:ncol(centroidMat)) {
    sampSub.1 <- which(clusterRes$subtypeCall == names(centroidMat)[j])
    if(length(sampSub.1)!=0) {
      datTmp <-  clusterRes[sampSub.1,]
      sampSrt <- rownames(datTmp)[order(datTmp$r)]
      sampSrtAll <- c(sampSrtAll,sampSrt)
    }
  }
  clusterRes$order <- match(sampSrtAll,rownames(clusterRes))
  dataSet$tmp$clusterRes <- clusterRes
  dataSet$tmp$subtypeCall <- clusterRes$subtypeCall
  
  if(FALSE){
  if(lowVarFlt){
    # filter rows with smallest variance
    tmpEx <- ex
    #tmpEx[tmpEx < 0.5] <- 0
    featVars <- apply(tmpEx, 1, sd)
    # idxFlt <- which(featVars==min(featVars))
    idxFlt <- which(featVars < quantile(featVars,0.2))
    ex <- ex[-idxFlt, ]
    dataSet$tmp$featInfo <- dataSet$tmp$featInfo[-idxFlt]
  }
  
  if(Rversion == "R-4.X.X") {
    clusterAlg="km"
  } else {
    clusterAlg="kmdist"
  }
  
  ### row
  dataSet$tmp$Rowv <- as.dendrogram(
    ConsensusClusterPlus::ConsensusClusterPlus(
      d = t(as.matrix(ex)),
      seed = 1234, pFeature = 0.8, pItem = 0.8,
      maxK = 5+1, reps=200, distance=distance,
      clusterAlg=clusterAlg)[[5]]$consensusTree)
  }
  
  # save subtype to sampInfo
  dataSet$sampInfo$tmpCluster <- NA
  dataSet$sampInfo$tmpCluster[sampSub] <- dataSet$tmp$subtypeCall
  names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- dataSet$tmp$mainLab
  
  # color samples
  ColSideColors <- c(subtypeColList[[i]])[match(dataSet$tmp$subtypeCall, subtypeList[[i]])]
  dataSet$tmp$ColSideColors <- ColSideColors
  
  # rename tmp 
  dataSet[[dataSet$tmp$mainLab]] <- dataSet$tmp
  dataSet$tmp <- NULL
  
  return(dataSet)
}