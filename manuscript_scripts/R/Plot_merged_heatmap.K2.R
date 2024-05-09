Plot_merged_heatmap.K2 <- function(rDataName, dataSet, cafSubtype, sampSub, sampSubLabel = "") {
  
  if(sampSubLabel == "") {
    pdf(paste("../heatmap/",rDataName,".K2.pdf",sep=""))
  } else {
    pdf(paste("../heatmap/",rDataName,".",sampSubLabel,".K2.pdf",sep=""))
  }
  
  
  ColSideColors0 <- data.frame(PurIST = Get_subtype_color(cafSubtype$Subtype$PurIST,"PurIST"),                     
                               PurIST.prob = colorRampPalette(c('snow2', 'black'))(length(cafSubtype$Subtype$PurIST.prob))[rank(cafSubtype$Subtype$PurIST.prob)],
                               stringsAsFactors = FALSE)
  
  ColSideColors <- cbind(ColSideColors0[sampSub,],
                         data.frame(MS.unscaled = cafSubtype$MS.unscaled$ColSideColors,
                                    Elyada_CAF.unscaled = cafSubtype$Elyada_CAF.unscaled$ColSideColors,
                                    SCISSORS_CAF_K2_top10.unscaled = cafSubtype$SCISSORS_CAF_K2_top10.unscaled$ColSideColors,
                                    stringsAsFactors = FALSE))
  
  ### purist
  TSPgenes <- subtypeGeneList[[1]]
  featSet <- match(TSPgenes$geneSymbol,dataSet$featInfo$SYMBOL)
  sampOrd <- order(cafSubtype$Subtype$PurIST.prob[sampSub], decreasing = TRUE)
  TSPgeneMat <- dataSet$ex[featSet, sampOrd]
  TSPgeneMat <- apply(TSPgeneMat, 2, rank)
  Plot_heatmap_ordered(TSPgeneMat, ColSideColors[sampOrd,], TSPgenes, "PurIST TSP")
  
  Plot_heatmap_CC(cafSubtype$MS.unscaled, "MS", ColSideColors, cafSubtype$MS.unscaled$mainLab.2)
  Plot_heatmap_CC(cafSubtype$Elyada_CAF.unscaled, "Elyada_CAF", ColSideColors, cafSubtype$Elyada_CAF.unscaled$mainLab.2)
  Plot_heatmap_CC(cafSubtype$SCISSORS_CAF_K2_top10.unscaled, "SCISSORS_CAF_K2_top10", ColSideColors, cafSubtype$SCISSORS_CAF_K2_top10.unscaled$mainLab.2)
  dev.off()
}