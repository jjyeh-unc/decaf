Plot_merged_heatmap.K4 <- function(rDataName, dataSet, cafSubtype, sampSub, sampSubLabel = "") {
  
  if(sampSubLabel == "") {
    pdf(paste("../heatmap/",rDataName,".K4.pdf",sep=""))
  } else {
    pdf(paste("../heatmap/",rDataName,".",sampSubLabel,".K4.pdf",sep=""))
  }
  
  
  ColSideColors0 <- data.frame(PurISS.71322 = Get_subtype_color(cafSubtype$Subtype$PurISS.71322,"PurISS.71322"),   
                               PurISS_graded.71322 = Get_subtype_color(cafSubtype$Subtype$PurISS_graded.71322,"PurISS_graded.71322"),   
                               PurISS.prob.71322 = colorRampPalette(c('snow2', 'black'))(length(cafSubtype$Subtype$PurISS.prob.71322))[rank(cafSubtype$Subtype$PurISS.prob.71322)],
                               stringsAsFactors = FALSE)
  
  ColSideColors <- cbind(ColSideColors0[sampSub,],
                         data.frame(SCISSORS_CAF_top25 = cafSubtype$SCISSORS_CAF_K2_top25.unscaled$ColSideColors,
                                    stringsAsFactors = FALSE))
  
  ### purist
  TSPgenes <- subtypeGeneList[[27]]
  featSet <- match(TSPgenes$geneSymbol,dataSet$featInfo$SYMBOL)
  TSPgeneMat <- dataSet$ex[featSet, sampSub]
  sampOrd <- order(cafSubtype$Subtype$PurISS.prob.71322[sampSub], decreasing = TRUE)
  TSPgeneMat <- TSPgeneMat[, sampOrd]
  TSPgeneMat <- apply(TSPgeneMat, 2, rank)
  Plot_heatmap_ordered(TSPgeneMat, ColSideColors[sampOrd,], TSPgenes, "PurISS.71322 TSP")
  
  Plot_heatmap_CC(cafSubtype$SCISSORS_CAF_K2_top25.unscaled, "SCISSORS_CAF_K2_top25", ColSideColors, cafSubtype$SCISSORS_CAF_K2_top25.unscaled$mainLab.2)
  dev.off()
}