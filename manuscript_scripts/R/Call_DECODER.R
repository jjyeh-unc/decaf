Call_DECODER <- function(dataSet) {
  if(length(dataSet$metadata$log.transformed) == 1){
    if(!(dataSet$metadata$log.transformed)){
      tmpEx <- log2(dataSet$ex+1)
    } else{
      tmpEx <- dataSet$ex
    }
  }else{
    tmpEx <- log2(dataSet$ex+1)
  }
  rownames(tmpEx) <- make.names(dataSet$featInfo$SYMBOL, unique = T)
  sampleWeight <- Decon_single_sample("TCGA_RNAseq_PAAD",
                                       tmpEx,
                                       "geneSymbol")
  #dataSet$decoderWeights  <- Norm_PDAC_weights(sampleWeights)
  #dataSet$sampInfo$DECODER.tumorCall <- NA
  #dataSet$sampInfo$DECODER.tumorCall[dataSet$decoderWeights$bcRatio >= 1] <- "Basal-like"
  #dataSet$sampInfo$DECODER.tumorCall[dataSet$decoderWeights$bcRatio < 1] <- "Classical"
  
  sampleWeight <- sampleWeight[,c(9,5,4,7,2,1,3)]
  sampleWeightNorm <- data.frame(sweep(sampleWeight, 1, rowSums(sampleWeight), FUN="/"))
  names(sampleWeightNorm) <- c("BasalTumor","ClassicalTumor",
                               "ActivatedStroma","NormalStroma",
                               "Immune","Endocrine","Exocrine")
  
  sampleWeightNorm$bcRatio <- (sampleWeightNorm$BasalTumor+0.01)/(sampleWeightNorm$ClassicalTumor+0.01)
  sampleWeightNorm$bcRatio[sampleWeightNorm$bcRatio > 5] <- 5
  sampleWeightNorm$bcDiff <- sampleWeightNorm$BasalTumor-sampleWeightNorm$ClassicalTumor
  sampleWeightNorm$bcSum <- sampleWeightNorm$BasalTumor+sampleWeightNorm$ClassicalTumor
  
  sampleWeightNorm$anRatio <- (sampleWeightNorm$ActivatedStroma+0.01)/(sampleWeightNorm$NormalStroma+0.01)
  sampleWeightNorm$anRatio[sampleWeightNorm$anRatio > 5] <- 5
  sampleWeightNorm$anDiff <- sampleWeightNorm$ActivatedStroma-sampleWeightNorm$NormalStroma
  sampleWeightNorm$anSum <- sampleWeightNorm$ActivatedStroma+sampleWeightNorm$NormalStroma
  
  sampleWeightNorm$aiRatio <- (sampleWeightNorm$ActivatedStroma+0.01)/(sampleWeightNorm$Immune+0.01)
  sampleWeightNorm$aiRatio[sampleWeightNorm$aiRatio > 5] <- 5
  sampleWeightNorm$aiDiff <- sampleWeightNorm$ActivatedStroma-sampleWeightNorm$Immune
  sampleWeightNorm$aiSum <- sampleWeightNorm$ActivatedStroma+sampleWeightNorm$Immune
  
  sampleWeightNorm$aniRatio <- (sampleWeightNorm$ActivatedStroma+0.01)/(sampleWeightNorm$NormalStroma+sampleWeightNorm$Immune+0.01)
  sampleWeightNorm$aniRatio[sampleWeightNorm$aniRatio > 5] <- 5
  sampleWeightNorm$aniDiff <- sampleWeightNorm$ActivatedStroma-sampleWeightNorm$NormalStroma-sampleWeightNorm$Immune
  sampleWeightNorm$aniSum <- sampleWeightNorm$ActivatedStroma+sampleWeightNorm$NormalStroma+sampleWeightNorm$Immune
  
  dataSet$decoderWeights  <- sampleWeightNorm
  
  #decoderWeights <- dataSet$decoderWeights 
  #names(decoderWeights) <- paste("DECODER.",names(decoderWeights),sep="")
  #dataSet$sampInfo <- cbind(dataSet$sampInfo, decoderWeights)
  return(dataSet)
}