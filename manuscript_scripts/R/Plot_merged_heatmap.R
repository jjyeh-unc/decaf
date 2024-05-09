Plot_merged_heatmap <- function(rDataName, dataSet, cafSubtype, sampSub, sampSubLabel = "") {
  
  if(sampSubLabel == "") {
    pdf(paste("../heatmap/",rDataName,".pdf",sep=""))
  } else {
    pdf(paste("../heatmap/",rDataName,".",sampSubLabel,".pdf",sep=""))
  }
  
  
  ColSideColors0 <- data.frame(PurIST = Get_subtype_color(cafSubtype$Subtype$PurIST,"PurIST"),                     
                               PurIST.prob = colorRampPalette(c('snow2', 'black'))(length(cafSubtype$Subtype$PurIST.prob))[rank(cafSubtype$Subtype$PurIST.prob)],
                               stringsAsFactors = FALSE)
  
  ColSideColors <- cbind(ColSideColors0[sampSub,],
                         data.frame(MS.unscaled = cafSubtype$MS.unscaled$ColSideColors,
                                    MS.scaled = cafSubtype$MS.scaled$ColSideColors,
                                    Elyada_CAF.unscaled = cafSubtype$Elyada_CAF.unscaled$ColSideColors,
                                    Elyada_CAF.scaled = cafSubtype$Elyada_CAF.scaled$ColSideColors,
                                    SCISSORS_CAF_K2_top25.unscaled = cafSubtype$SCISSORS_CAF_K2_top25.unscaled$ColSideColors,
                                    SCISSORS_CAF_K2_top25.scaled = cafSubtype$SCISSORS_CAF_K2_top25.scaled$ColSideColors,
                                    SCISSORS_CAF_K2_top10.unscaled = cafSubtype$SCISSORS_CAF_K2_top10.unscaled$ColSideColors,
                                    SCISSORS_CAF_K2_top10.scaled = cafSubtype$SCISSORS_CAF_K2_top10.scaled$ColSideColors,
                                    SCISSORS_CAF_top25.unscaled = cafSubtype$SCISSORS_CAF_top25.unscaled$ColSideColors,
                                    SCISSORS_CAF_top10.unscaled = cafSubtype$SCISSORS_CAF_top10.unscaled$ColSideColors,
                                    SCISSORS_panCAF_vs_peri_top25.unscaled = cafSubtype$SCISSORS_panCAF_vs_peri_top25.unscaled$ColSideColors,
                                    SCISSORS_panCAF_vs_peri_top10.unscaled = cafSubtype$SCISSORS_panCAF_vs_peri_top10.unscaled$ColSideColors,
                                    stringsAsFactors = FALSE))
  
  ### purist
  TSPgenes <- subtypeGeneList[[1]]
  featSet <- match(TSPgenes$geneSymbol,dataSet$featInfo$SYMBOL)
  sampOrd <- order(cafSubtype$Subtype$PurIST.prob[sampSub], decreasing = TRUE)
  TSPgeneMat <- dataSet$ex[featSet, sampOrd]
  TSPgeneMat <- apply(TSPgeneMat, 2, rank)
  Plot_heatmap_ordered(TSPgeneMat, ColSideColors[sampOrd,], TSPgenes, "PurIST TSP")
  
  Plot_heatmap_CC(cafSubtype$MS.unscaled, "MS", ColSideColors, cafSubtype$MS.unscaled$mainLab.2)
  Plot_heatmap_CC(cafSubtype$MS.scaled, "MS", ColSideColors, cafSubtype$MS.scaled$mainLab.2)
  Plot_heatmap_CC(cafSubtype$Elyada_CAF.unscaled, "Elyada_CAF", ColSideColors, cafSubtype$Elyada_CAF.unscaled$mainLab.2)
  Plot_heatmap_CC(cafSubtype$Elyada_CAF.scaled, "Elyada_CAF", ColSideColors, cafSubtype$Elyada_CAF.scaled$mainLab.2)
  Plot_heatmap_CC(cafSubtype$SCISSORS_CAF_K2_top25.unscaled, "SCISSORS_CAF_K2_top25", ColSideColors, cafSubtype$SCISSORS_CAF_K2_top25.unscaled$mainLab.2)
  Plot_heatmap_CC(cafSubtype$SCISSORS_CAF_K2_top25.scaled, "SCISSORS_CAF_K2_top25", ColSideColors, cafSubtype$SCISSORS_CAF_K2_top25.scaled$mainLab.2)
  Plot_heatmap_CC(cafSubtype$SCISSORS_CAF_K2_top10.unscaled, "SCISSORS_CAF_K2_top10", ColSideColors, cafSubtype$SCISSORS_CAF_K2_top10.unscaled$mainLab.2)
  Plot_heatmap_CC(cafSubtype$SCISSORS_CAF_K2_top10.scaled, "SCISSORS_CAF_K2_top10", ColSideColors, cafSubtype$SCISSORS_CAF_K2_top10.scaled$mainLab.2)
  Plot_heatmap_CC(cafSubtype$SCISSORS_CAF_top25.unscaled, "SCISSORS_CAF_top25", ColSideColors, cafSubtype$SCISSORS_CAF_top25.unscaled$mainLab.2)
  Plot_heatmap_CC(cafSubtype$SCISSORS_CAF_top10.unscaled, "SCISSORS_CAF_top10", ColSideColors, cafSubtype$SCISSORS_CAF_top10.unscaled$mainLab.2)
  Plot_heatmap_CC(cafSubtype$SCISSORS_panCAF_vs_peri_top25.unscaled, "SCISSORS_panCAF_vs_peri_top25", ColSideColors, cafSubtype$SCISSORS_panCAF_vs_peri_top25.unscaled$mainLab.2)
  Plot_heatmap_CC(cafSubtype$SCISSORS_panCAF_vs_peri_top10.unscaled, "SCISSORS_panCAF_vs_peri_top10", ColSideColors, cafSubtype$SCISSORS_panCAF_vs_peri_top10.unscaled$mainLab.2)
  dev.off()
}