save_subtype <- function(dataSet, tmpCalls, schemaName, sampSub){
  if(length(dataSet$subtypeCalls) == 0){
    dataSet$subtypeCalls <- data.frame(tmp = matrix(nrow = (nrow(dataSet$sampInfo))),
                                       stringsAsFactors = FALSE)
    colnames(dataSet$subtypeCalls) <- schemaName
    dataSet$subtypeCalls[sampSub, schemaName] <- tmpCalls
  } else{
    dataSet$subtypeCalls$schemaName <- rep_len(NA, length.out = nrow(dataSet$sampInfo))
    names(dataSet$subtypeCalls)[which(names(dataSet$subtypeCalls) == "schemaName")] <- schemaName
    dataSet$subtypeCalls[sampSub, schemaName] <- tmpCalls
  }
  return(dataSet)
}

