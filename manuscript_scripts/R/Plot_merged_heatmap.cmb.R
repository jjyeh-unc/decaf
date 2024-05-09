Plot_merged_heatmap.cmb <- function(rDataName, dataSet, cafSubtype, cafSubtype.July, cafSubtype.July_pi, sampSub) {
  
  pdf(paste("../heatmap/",rDataName,".cmb.pdf",sep=""))

  ColSideColors0 <- data.frame(PurIST = Get_subtype_color(cafSubtype$Subtype$PurIST,"PurIST"),   
                               PurIST_graded = Get_subtype_color(cafSubtype$Subtype$PurIST_graded,"PurIST_graded"),   
                               PurIST.prob = colorRampPalette(c('snow2', 'black'))(length(cafSubtype$Subtype$PurIST.prob))[rank(cafSubtype$Subtype$PurIST.prob)],
                               PurISS.71322 = Get_subtype_color(cafSubtype$Subtype$PurISS.71322,"PurISS.71322"),   
                               PurISS_graded.71322 = Get_subtype_color(cafSubtype$Subtype$PurISS_graded.71322,"PurISS_graded.71322"),   
                               PurISS.prob.71322 = colorRampPalette(c('snow2', 'black'))(length(cafSubtype$Subtype$PurISS.prob.71322))[rank(cafSubtype$Subtype$PurISS.prob.71322)],
                               stringsAsFactors = FALSE)
  if(setequal(sampSub, cafSubtype.July$MS.unscaled$sampSub)) {
    ColSideColors <- cbind(ColSideColors0[sampSub,],
                           data.frame(SCISSORS_top25_Nov = cafSubtype$SCISSORS_CAF_K2_top25.unscaled$ColSideColors,
                                      SCISSORS_top25_July_pi = cafSubtype.July_pi$SCISSORS_CAF_K2_top25.unscaled$ColSideColors,
                                      SCISSORS_top10_July_pi = cafSubtype.July_pi$SCISSORS_CAF_K2_top10.unscaled$ColSideColors,
                                      SCISSORS_top25_July = cafSubtype.July$SCISSORS_CAF_K2_top25.unscaled$ColSideColors,
                                      SCISSORS_top10_July = cafSubtype.July$SCISSORS_CAF_K2_top10.unscaled$ColSideColors,
                                      Elyada = cafSubtype.July$Elyada_CAF.unscaled$ColSideColors,
                                      MS = cafSubtype.July$MS.unscaled$ColSideColors,
                                      Puleo = cafSubtype$Puleo$ColSideColors,
                                      stringsAsFactors = FALSE))
    ColSideColors.1 <- ColSideColors
  } else {
    ColSideColors <- cbind(ColSideColors0[sampSub,],
                           data.frame(SCISSORS_top25_Nov = cafSubtype$SCISSORS_CAF_K2_top25.unscaled$ColSideColors,
                                      SCISSORS_top25_July_pi = cafSubtype.July_pi$SCISSORS_CAF_K2_top25.unscaled$ColSideColors,
                                      SCISSORS_top10_July_pi = cafSubtype.July_pi$SCISSORS_CAF_K2_top10.unscaled$ColSideColors,
                                      SCISSORS_top25_July = cafSubtype.July$SCISSORS_CAF_K2_top25.unscaled$ColSideColors[sampSub],
                                      SCISSORS_top10_July = cafSubtype.July$SCISSORS_CAF_K2_top10.unscaled$ColSideColors[sampSub],
                                      Elyada = cafSubtype.July$Elyada_CAF.unscaled$ColSideColors[sampSub],
                                      MS = cafSubtype.July$MS.unscaled$ColSideColors[sampSub],
                                      Puleo = cafSubtype$Puleo$ColSideColors,
                                      stringsAsFactors = FALSE))
    
    ColSideColors.1 <- cbind(ColSideColors0,
                             data.frame(#SCISSORS_top25_Nov = cafSubtype$SCISSORS_CAF_K2_top25.unscaled$ColSideColors,
                                        SCISSORS_top25_July = cafSubtype.July$SCISSORS_CAF_K2_top25.unscaled$ColSideColors,
                                        SCISSORS_top10_July = cafSubtype.July$SCISSORS_CAF_K2_top10.unscaled$ColSideColors,
                                        Elyada = cafSubtype.July$Elyada_CAF.unscaled$ColSideColors,
                                        MS = cafSubtype.July$MS.unscaled$ColSideColors,
                                        #Puleo = cafSubtype$Puleo$ColSideColors,
                                        stringsAsFactors = FALSE))
  }
  
  ### puriss
  TSPgenes <- subtypeGeneList[[27]]
  featSet <- match(TSPgenes$geneSymbol,dataSet$featInfo$SYMBOL)
  TSPgeneMat <- dataSet$ex[featSet, sampSub]
  sampOrd <- order(cafSubtype$Subtype$PurISS.prob.71322[sampSub], decreasing = TRUE)
  TSPgeneMat <- TSPgeneMat[, sampOrd]
  TSPgeneMat <- apply(TSPgeneMat, 2, rank)
  Plot_heatmap_ordered(TSPgeneMat, ColSideColors[sampOrd,], TSPgenes, "PurISS.71322 TSP")
  
  Plot_heatmap_CC(cafSubtype$SCISSORS_CAF_K2_top25.unscaled, "SCISSORS_CAF_K2_top25", ColSideColors, paste("Nov\n",cafSubtype$SCISSORS_CAF_K2_top25.unscaled$mainLab.2,sep=""))
  Plot_heatmap_CC(cafSubtype.July_pi$SCISSORS_CAF_K2_top25.unscaled, "SCISSORS_CAF_K2_top25", ColSideColors, paste("July_pi\n",cafSubtype.July_pi$SCISSORS_CAF_K2_top25.unscaled$mainLab.2,sep=""))
  Plot_heatmap_CC(cafSubtype.July_pi$SCISSORS_CAF_K2_top10.unscaled, "SCISSORS_CAF_K2_top10", ColSideColors, paste("July_pi\n",cafSubtype.July_pi$SCISSORS_CAF_K2_top10.unscaled$mainLab.2,sep=""))
  Plot_heatmap_CC(cafSubtype.July$SCISSORS_CAF_K2_top25.unscaled, "SCISSORS_CAF_K2_top25", ColSideColors.1, paste("July\n",cafSubtype.July$SCISSORS_CAF_K2_top25.unscaled$mainLab.2,sep=""))
  Plot_heatmap_CC(cafSubtype.July$SCISSORS_CAF_K2_top10.unscaled, "SCISSORS_CAF_K2_top10", ColSideColors.1, paste("July\n",cafSubtype.July$SCISSORS_CAF_K2_top10.unscaled$mainLab.2,sep=""))
  Plot_heatmap_CC(cafSubtype.July$Elyada_CAF.unscaled, "Elyada_CAF", ColSideColors.1, paste("July\n",cafSubtype.July$Elyada_CAF.unscaled$mainLab.2,sep=""))
  Plot_heatmap_CC(cafSubtype.July$MS.unscaled, "MS", ColSideColors.1, paste("July\n",cafSubtype.July$MS.unscaled$mainLab.2,sep=""))
  dev.off()
}