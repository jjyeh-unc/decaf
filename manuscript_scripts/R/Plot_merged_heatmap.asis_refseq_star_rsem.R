Plot_merged_heatmap.asis_refseq_star_rsem <- function(rDataName, dataSet, cafSubtype, sampSub) {
  
  pdf(paste("../heatmap/",rDataName,".asis_refseq_star_rsem.pdf",sep=""))

  ColSideColors0 <- data.frame(PurIST = Get_subtype_color(cafSubtype$Subtype$PurIST,"PurIST"),   
                               PurIST_graded = Get_subtype_color(cafSubtype$Subtype$PurIST_graded,"PurIST_graded"),   
                               PurIST.prob = colorRampPalette(c('snow2', 'black'))(length(cafSubtype$Subtype$PurIST.prob))[rank(cafSubtype$Subtype$PurIST.prob)],
                               stringsAsFactors = FALSE)
  ColSideColors <- cbind(ColSideColors0[sampSub,],
                         data.frame(SCISSORS_vK = c("coral","peachpuff1","grey","palegreen3","darkgreen","black")[match(cafSubtype$Subtype$SCISSORS_CAF_K2_top25.vK , c("iCAF","Mixed.iCAF","Mixed","Mixed.myCAF","myCAF","Absent"))[sampSub]],
                                    SCISSORS_K2 = cafSubtype$SCISSORS_CAF_K2_top25.unscaled.K2$ColSideColors,
                                    Elyada = cafSubtype$Elyada_CAF.unscaled.K2$ColSideColors,
                                    MS = cafSubtype$MS.unscaled.K2$ColSideColors,
                                    Maurer = cafSubtype$Maurer.unscaled.K2$ColSideColors,
                                    Puleo = cafSubtype$Puleo$ColSideColors,
                                    stringsAsFactors = FALSE))
  ### purist
  TSPgenes <- subtypeGeneList[[1]]
  featSet <- match(TSPgenes$geneSymbol,dataSet$featInfo$SYMBOL)
  sampOrd <- order(cafSubtype$Subtype$PurIST.prob[sampSub], decreasing = TRUE)
  TSPgeneMat <- dataSet$ex[featSet, sampOrd]
  TSPgeneMat <- apply(TSPgeneMat, 2, rank)
  Plot_heatmap_ordered(TSPgeneMat, ColSideColors[sampOrd,], TSPgenes, "PurIST TSP")
  
  for (tmp in names(cafSubtype)) {
    if(tmp %in% c("Subtype","Puleo")) {next}
    Plot_heatmap_CC(cafSubtype[[tmp]], cafSubtype[[tmp]]$schema, ColSideColors, cafSubtype[[tmp]]$mainLab.2)
  } 
  dev.off()
}