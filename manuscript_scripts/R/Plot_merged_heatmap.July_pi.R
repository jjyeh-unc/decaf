Plot_merged_heatmap.July_pi <- function(rDataName, dataSet, cafSubtype, sampSub) {
  
  pdf(paste("../heatmap/",rDataName,".July_pi.pdf",sep=""))
  
  ColSideColors0 <- data.frame(PurIST = Get_subtype_color(cafSubtype$Subtype$PurIST,"PurIST"),   
                               PurIST_graded = Get_subtype_color(cafSubtype$Subtype$PurIST_graded,"PurIST_graded"),   
                               PurIST.prob = colorRampPalette(c('snow2', 'black'))(length(cafSubtype$Subtype$PurIST.prob))[rank(cafSubtype$Subtype$PurIST.prob)],
                               PurISS.71322 = Get_subtype_color(cafSubtype$Subtype$PurISS.71322,"PurISS.71322"),   
                               PurISS_graded.71322 = Get_subtype_color(cafSubtype$Subtype$PurISS_graded.71322,"PurISS_graded.71322"),   
                               PurISS.prob.71322 = colorRampPalette(c('snow2', 'black'))(length(cafSubtype$Subtype$PurISS.prob.71322))[rank(cafSubtype$Subtype$PurISS.prob.71322)],
                               stringsAsFactors = FALSE)
  
  ColSideColors <- cbind(ColSideColors0[sampSub,],
                         data.frame(SCISSORS_CAF_top25 = cafSubtype$SCISSORS_CAF_K2_top25.unscaled$ColSideColors,
                                    SCISSORS_CAF_top10 = cafSubtype$SCISSORS_CAF_K2_top10.unscaled$ColSideColors,
                                    stringsAsFactors = FALSE))
  
  ### puriss
  TSPgenes <- subtypeGeneList[[27]]
  featSet <- match(TSPgenes$geneSymbol,dataSet$featInfo$SYMBOL)
  TSPgeneMat <- dataSet$ex[featSet, sampSub]
  sampOrd <- order(cafSubtype$Subtype$PurISS.prob.71322[sampSub], decreasing = TRUE)
  TSPgeneMat <- TSPgeneMat[, sampOrd]
  TSPgeneMat <- apply(TSPgeneMat, 2, rank)
  Plot_heatmap_ordered(TSPgeneMat, ColSideColors[sampOrd,], TSPgenes, "PurISS.71322 TSP")
  
  Plot_heatmap_CC(cafSubtype$SCISSORS_CAF_K2_top25.unscaled, "SCISSORS_CAF_K2_top25", ColSideColors, cafSubtype$SCISSORS_CAF_K2_top25.unscaled$mainLab.2)
  Plot_heatmap_CC(cafSubtype$SCISSORS_CAF_K2_top10.unscaled, "SCISSORS_CAF_K2_top10", ColSideColors, cafSubtype$SCISSORS_CAF_K2_top10.unscaled$mainLab.2)
  dev.off()
}