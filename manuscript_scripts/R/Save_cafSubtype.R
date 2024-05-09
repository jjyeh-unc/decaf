Save_cafSubtype <- function(cafSubtype, dataSet, tmpSchema, tmpK, tmpRowScale, tmpVersion="") {
  
  tmpRowScaleLab <- c("unscaled", "scaled")[(tmpRowScale + 1)]
  tmpMainLab.2 <- paste(tmpSchema,".",tmpRowScaleLab,".","K",tmpK,sep="")
  tmpMainLab <- paste(tmpSchema,".",tmpRowScaleLab,sep="")
  
  if(tmpVersion != "") {
    tmpMainLab.0 <- paste(tmpSchema,".",tmpVersion,sep="")
  } else {
    tmpMainLab.0 <- tmpSchema
  }
  
  
  cafSubtype[[tmpMainLab.2]] <- dataSet[[tmpMainLab.2]]
  cafSubtype$Subtype[tmpMainLab.0] <- dataSet$sampInfo[[tmpMainLab.2]]
  cafSubtype$Subtype[paste(tmpMainLab.0,".","classifier",sep="")] <- FALSE
  cafSubtype$Subtype[which(!is.na(cafSubtype$Subtype[tmpMainLab.0] )),paste(tmpMainLab.0,".","classifier",sep="")] <- TRUE
  return(cafSubtype)
}