Call_consensus <- function(dataSet, schemaTmp, sampSub, 
                           rowK = 2, colK = 2, 
                           rowScale = FALSE, colScale = FALSE, lowVarFlt = 0,
                           clusterAlg = "km", distance = "euclidean", Rversion = "R-4.X.X") {
  
  dataSet$tmp <- list()
  tmpRowScaleLab <- c("unscaled", "scaled")[(rowScale + 1)]
  dataSet$tmp$mainLab.2 <- paste(schemaTmp,".",tmpRowScaleLab,".","K",colK,sep="")
  dataSet$tmp$mainLab <- paste(schemaTmp,".",tmpRowScaleLab,sep="")
  dataSet$tmp$schema <- schemaTmp
  dataSet$tmp$sampSub <- sampSub
  dataSet$tmp$rowK <- rowK
  dataSet$tmp$colK <- colK
  dataSet$tmp$rowScale <- rowScale
  dataSet$tmp$colScale <- colScale
  dataSet$tmp$clusterAlg <- clusterAlg
  dataSet$tmp$distance <- distance
  dataSet$tmp$Rversion <- Rversion
  
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
    datCol <- as.dist(1-cor( as.matrix(dataSet$tmp$ex),method="pearson"))
    datRow <- as.dist(1-cor(t(as.matrix(dataSet$tmp$ex)),method="pearson"))
  } else {
    datCol <- as.matrix(dataSet$tmp$ex)
    datRow <- t(as.matrix(dataSet$tmp$ex))
  }
  
  dataSet$tmp$Colv <- as.dendrogram(
    ConsensusClusterPlus::ConsensusClusterPlus(
      d = datCol,
      seed = 9999, pFeature = 1, pItem = 0.8,
      maxK = colK+1, reps=1000, distance=distance,
      clusterAlg=clusterAlg)[[colK]]$consensusTree)
  dataSet$tmp$subtypeCall <- cutree(as.hclust(dataSet$tmp$Colv), k=colK)
  
  dataSet$tmp$Rowv <- as.dendrogram(
    ConsensusClusterPlus::ConsensusClusterPlus(
      d = datRow,
      seed = 9999, pFeature = 1, pItem = 0.8,
      maxK = rowK+1, reps=200, distance=distance,
      clusterAlg=clusterAlg)[[rowK]]$consensusTree)
  
  dataSet$tmp$sampID <- colnames(dataSet$ex)[sampSub]
  
  # plot heatmap -----------------------------------------------------------------------------------------------------
  tmp <- dataSet$tmp
  i <- which(schemaList %in% schemaTmp)
  genesTmp <- subtypeGeneList[[i]]
  idxGene <- match(tmp$featInfo,genesTmp$geneSymbol)
  geneInfo <- genesTmp[idxGene,]
  
  #tmp$subtype <- dataSet$sampInfo$tmpCluster[sampSub]
  tmp$ColSideColors <- c(brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 8, name = "Dark2"))[tmp$subtypeCall]
  
  if(ncol(as.matrix(tmp$ColSideColors))<=10) {
    ColSideColorsSize <- ncol(as.matrix(tmp$ColSideColors))
  } else {
    ColSideColorsSize <- 10
  }

   heatmap.3(main= dataSet$tmp$mainLab.2
            ,x = as.matrix(tmp$ex)
            ,Rowv = dataSet$tmp$Rowv
            ,Colv = dataSet$tmp$Colv
            ,dendrogram = "none" ,trace="none" ,scale="row"
            ,labRow = tmp$featInfo
            ,labCol = tmp$sampID
            ,col = colorRampPalette(c("darkblue","blue","white","red","darkred"))(99)
            ,ColSideColors = as.matrix(tmp$ColSideColors)
            ,ColSideColorsSize = ColSideColorsSize
            ,RowSideColors = t(as.matrix(geneInfo$Color))
            ,RowSideColorsSize=1
            ,lwid = c(1,5)
            ,lhei = c(1,3)
            ,margins =c(8,18)
            ,keysize=0.5
            ,cex.lab=0.4
            ,cexRow = 0.5
            ,cexCol = 0.6
  )
  
  legend(xy.coords(x=.8,y=.98),
         legend=c("Gene subtype", subtypeList[[i]],
                  "Cluster",as.character(1:20)),
         fill=c("white",
                subtypeColList[[i]],
                "white",
                brewer.pal(n = 12, name = "Set3"),
                brewer.pal(n = 8, name = "Dark2")),
         border=FALSE, bty="n",
         y.intersp = 1 , cex = 0.8)
  
  # rename tmp -----------------------------------------------------------------------------------------------------
  dataSet[[tmp$mainLab.2]] <- dataSet$tmp
  dataSet$tmp <- NULL
  
  return(dataSet)
}
