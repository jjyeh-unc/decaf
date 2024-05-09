Relabel_subtype <- function(dataSet, tmpSubtypeLab, tmpSchema, tmpK, tmpRowScale) {
  
  tmpRowScaleLab <- c("unscaled", "scaled")[(tmpRowScale + 1)]
  tmpMainLab.2 <- paste(tmpSchema,".",tmpRowScaleLab,".","K",tmpK,sep="")
  tmpMainLab <- paste(tmpSchema,".",tmpRowScaleLab,sep="")
  
  tmp <- dataSet[[tmpMainLab.2]]
  sampSub <- tmp$sampSub
  
  # find schema #
  i <- which(schemaList %in% tmp$schema)
  
  # assign subtype label
  tmp$subtypeCall.cluster <- tmp$subtypeCall
  tmp$subtypeCall <- tmpSubtypeLab[tmp$subtypeCall]
  
  # save subtype to sampInfo
  dataSet$sampInfo$tmpCluster <- NA
  dataSet$sampInfo$tmpCluster[sampSub] <- tmp$subtypeCall
  names(dataSet$sampInfo)[which(names(dataSet$sampInfo) %in% "tmpCluster")] <- tmp$mainLab.2
  
  # color samples
  ColSideColors <- c(subtypeColList[[i]])[match(tmp$subtypeCall, subtypeList[[i]])]
  ColSideColors[which(!(ColSideColors %in% c(subtypeColList[[i]])))] <- NA
  ColSideColors[which(tmp$subtypeCall %in% "Absent")] <- "black"
  ColSideColors[which(tmp$subtypeCall %in% "Mixed")] <- "grey"
  ColSideColors[which(tmp$subtypeCall %in% "Mixed.myCAF")] <- "#7FB17F"
  ColSideColors[which(tmp$subtypeCall %in% "Mixed.iCAF")] <- "#FFBFA7"
  tmp$ColSideColors <- ColSideColors
  tmp$subtypeCall.cluster = factor(tmp$subtypeCall.cluster,levels = c("iCAF", "Mixed.iCAF", "Mixed", "Mixed.myCAF", "myCAF"))
  
  
  
  cc_subtype_names = c("iCAF", "Mixed.iCAF","Mixed","Mixed.myCAF", "myCAF")
  
  
  tmp$subtypeCall = tmp$subtypeCall.cluster
  # plot heatmap
  Plot_heatmap_CC(tmp, tmpSchema, ColSideColors, tmp$mainLab.2)
  
  # assign tmp
  dataSet[[tmpMainLab.2]] <- tmp
  
  return(dataSet)
}
