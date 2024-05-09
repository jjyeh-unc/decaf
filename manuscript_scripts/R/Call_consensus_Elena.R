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
  
  
  
  colnames(tmp$ex) = paste("C",1:ncol(tmp$ex), sep = "_")
  rownames(tmp$ex) = tmp$featInfo
  
  # create annotation vector for the samples
  annotation_col = data.frame('Scissors CAF' =tmp$subtypeCall, "PurIST" = Moffitt_GEO_array.caf_subtype$Subtype$PurIST[Moffitt_GEO_array.caf_subtype$Subtype$SCISSORS_CAF_K2_top25.vK.classifier],
                              "Elyada CAF" = Moffitt_GEO_array.caf_subtype$Subtype$Elyada_CAF[Moffitt_GEO_array.caf_subtype$Subtype$SCISSORS_CAF_K2_top25.vK.classifier],
                              "Moffitt Stroma" = Moffitt_GEO_array.caf_subtype$Subtype$MS[Moffitt_GEO_array.caf_subtype$Subtype$SCISSORS_CAF_K2_top25.vK.classifier],check.names = F)
  rownames(annotation_col) = colnames(tmp$ex)
  
  annotation_row = data.frame('Up Subtype' =rep("iCAF", length(geneInfo$myCAF)), check.names = F)
  annotation_row$`Up Subtype`[geneInfo$myCAF] = "myCAF"
  rownames(annotation_row) = rownames(tmp$ex)
  
  
  colscale2 = colorRampPalette(c("coral","white", "darkgreen"))(5)
  
  # create annotation color vector
  index_caf = which(schemaList == "SCISSORS_CAF_K2_top10")
  index_PurIST = which(schemaList == "PurIST")
  index_ms = which(schemaList == "MS_K2")
  annotation_col_colors = list(
    'Scissors CAF' = colscale2,
    "Up Subtype" = subtypeColList[[index_caf]],
    "PurIST" = subtypeColList[[index_PurIST]],
    "Elyada CAF" = subtypeColList[[index_caf]],
    "Moffitt Stroma" = subtypeColList[[index_ms]]
    
  )
  
  cc_subtype_names = c("iCAF", "Mixed.iCAF","Mixed","Mixed.myCAF", "myCAF")
  
  
  names(annotation_col_colors$'Scissors CAF')  = factor(cc_subtype_names, levels = c("iCAF", "Mixed.iCAF","Mixed","Mixed.myCAF", "myCAF"))
  names(annotation_col_colors$`Up Subtype`) =  subtypeList[[index_caf]]
  names(annotation_col_colors$"PurIST") =  subtypeList[[index_PurIST]]
  names(annotation_col_colors$"Elyada CAF") =  subtypeList[[index_caf]]
  names(annotation_col_colors$"Moffitt Stroma" ) =  subtypeList[[index_ms]]
  
  
  
  
  ## Sort columns and rows of heatmap
  sort = order(PurISS_data$Pred_prob_myCAF[samples_subset])
  tsp = match(as.vector(t(classifier$TSPs[classifier$fit$beta[-1]!=0,])),rownames(TSPgeneMat))
  
  ## Classifier gene annotation
  annotation_row = data.frame('Up Subtype' = factor(rev(rep(subtypeList[[index_caf]]   ,length(tsp)/2))),check.names = F)
  rownames(annotation_row) = rownames(TSPgeneMat)[tsp]
  
  
  ## Generate heatmap
  heatmap = pheatmap(TSPgeneMat[tsp,sort], annotation_col = annotation_col[sort,],annotation_row = annotation_row, 
                     annotation_colors = annotation_col_colors,show_colnames = F, cluster_cols = F,cluster_rows = F,
                     annotation_names_row = F,width = 8, height = 8)
  
  
  
  ColSideColors_ann = data.frame("Scissors CAF" = tmp$subtypeCall)
  row.names(ColSideColors_ann) = colnames(tmp$ex)
  
  z = q$colInd
  pheatmap(tmp$ex[q$rowInd,z[c(47:94,24:46, 1:23, 136:145, 95:135)]], scale = "row", annotation_col = annotation_col, annotation_row = annotation_row,
           color = colorRampPalette(c("darkblue","blue","white","red","darkred"))(99),
           annotation_colors = annotation_col_colors,
           show_colnames = F, cluster_cols = F,cluster_rows = F, annotation_names_row = F, 
           main = "Moffitt GEO array")
  
  pdf("heatmap_Scissors_CAF_PurIST_Elyada_MS_Moffitt_heatmap.pdf")
  pheatmap(tmp$ex[q$rowInd,z[c(47:94,24:46, 1:23, 136:145, 95:135)]], scale = "row", annotation_col = annotation_col, annotation_row = annotation_row,
           color = colorRampPalette(c("darkblue","blue","white","red","darkred"))(99),
           annotation_colors = annotation_col_colors,
           show_colnames = F, cluster_cols = F,cluster_rows = F, annotation_names_row = F, 
           main = "Moffitt GEO array")
  dev.off()
  pheatmap(tmp$ex[q$rowInd,], annotation_col = annotation_col,annotation_row = annotation_row, 
           annotation_colors = annotation_col_colors,show_colnames = F, cluster_cols = F,cluster_rows = F,
           annotation_names_row = F,width = 8, height = 8, scale = "row")
  
  
  
  
  pheatmap(tmp$ex)
  tmp$ex
  
   heatmap.3(main= tmp$mainLab.2
            ,x = as.matrix(tmp$ex)[41:1,]
            ,dendrogram = "none" ,trace="none" ,scale="row"
            ,labRow = tmp$featInfo
            ,labCol = tmp$sampID
            ,col = colorRampPalette(c("darkblue","blue","white","red","darkred"))(99)
            ,ColSideColors = as.matrix(tmp$ColSideColors)
            ,ColSideColorsSize = ColSideColorsSize
            ,RowSideColors = t(as.matrix(geneInfo$Color[41:1]))
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
