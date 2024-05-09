Plot_heatmap_CC_CAF <- function(tmp = cafSubtype[[tmp]], schemaTmp, ColSideColors, mainLab=schemaTmp) {
  i <- which(schemaList %in% schemaTmp)
  genesTmp <- subtypeGeneList[[i]]
  idxGene <- match(tmp$featInfo,genesTmp$geneSymbol)
  geneInfo <- genesTmp[idxGene,]
  
  if(ncol(as.matrix(ColSideColors))<=15) {
    ColSideColorsSize <- ncol(as.matrix(ColSideColors))
  } else {
    ColSideColorsSize <- 15
  }
  
  heatmap.3(main= mainLab
            ,x = as.matrix(tmp$ex)
            ,Rowv = tmp$Rowv
            ,Colv = tmp$Colv
            ,dendrogram = "both" ,trace="none" ,scale="row"
            ,labRow = tmp$featInfo
            ,labCol = ""
            ,col = colorRampPalette(c("darkblue","blue","white","red","darkred"))(99)
            ,ColSideColors = as.matrix(ColSideColors)
            ,ColSideColorsSize = ColSideColorsSize
            ,RowSideColors = t(as.matrix(geneInfo$Color))
            ,RowSideColorsSize=1
            ,lwid = c(1,5)
            ,lhei = c(1,3)
            ,margins =c(4,8)
            ,keysize=0.5
            ,cex.lab=0.4
            ,cexRow = 0.5
            ,cexCol = 0.6
  )
}
