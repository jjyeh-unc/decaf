Cluster_genes <- function(dataSet, schemaTmp, sampSub, 
                           rowK = 2, colK = 2, 
                           rowScale = FALSE, colScale = FALSE, lowVarFlt = 0,
                           clusterAlg = "km", distance = "euclidean", Rversion = "R-4.X.X") {
  dataSet$tmp <- list()
  if(length(dataSet$metadata$log.transformed) == 1){
    if(!(dataSet$metadata$log.transformed)){
      dataSet$tmp$ex <- log2(dataSet$ex+1)
    } else{
      dataSet$tmp$ex <- dataSet$ex
    }
  }else{
    dataSet$tmp$ex <- log2(dataSet$ex+1)
  }
  
  if(colScale){
    dataSet$tmp$ex <- preprocessCore::normalize.quantiles(as.matrix(dataSet$tmp$ex)) 
  }
  
  i <- which(schemaList %in% schemaTmp)
  geneList <- subtypeGeneList[[i]]
  geneList <- geneList$geneSymbol
  featSet <- which(dataSet$featInfo$SYMBOL %in% geneList)
  dataSet$tmp$ex <- dataSet$tmp$ex[featSet, sampSub]
  dataSet$tmp$featInfo <- dataSet$featInfo$SYMBOL[featSet]
  
  if(lowVarFlt){
    featVars <- apply(dataSet$tmp$ex, 1, sd)
    if(lowVarFlt == "TRUE") { # filter rows with 0 variance
      dataSet$tmp$ex <- dataSet$tmp$ex[featVars !=0, ]
      dataSet$tmp$featInfo <- dataSet$tmp$featInfo[featVars !=0]
    } else { # filter rows with smallest variance
      idxFlt <- which(featVars < quantile(featVars,lowVarFlt))
      dataSet$tmp$ex <- dataSet$tmp$ex[-idxFlt, ]
      dataSet$tmp$featInfo <- dataSet$tmp$featInfo[-idxFlt]
      }
  }
  
  if(rowScale){
    dataSet$tmp$ex <- t(scale(t(as.matrix(dataSet$tmp$ex))))
  }
  
  if( (Rversion == "R-4.X.X") & (clusterAlg=="km") ) {
    datRow <- as.dist(1-cor(t(as.matrix(dataSet$tmp$ex)),method="pearson"))
  } else {
    datRow <- t(as.matrix(dataSet$tmp$ex))
  }
  
  dataSet$tmp$Rowv <- as.dendrogram(
    ConsensusClusterPlus::ConsensusClusterPlus(
      d = datRow,
      seed = 9999, pFeature = 1, pItem = 0.8,
      maxK = rowK+1, reps=200, distance=distance,
      clusterAlg=clusterAlg)[[rowK]]$consensusTree)
  dataSet$tmp$featInfo <- cutree(as.hclust(dataSet$tmp$Rowv), k=rowK)
  
  dataSet$tmp$sampID <- dataSet$sampInfo[sampSub,1]
  return(dataSet)
}
