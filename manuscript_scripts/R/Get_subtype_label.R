Get_subtype_label <- function(callsTmp, schemaTmp) {
  i <- which(schemaList %in% schemaTmp)
  labTmp <- subtypeList[[i]]
  if(class(callsTmp) == "numeric") {
    labSamp <- labTmp[callsTmp]
  } else {
    labSamp <- labTmp[match(callsTmp,subtypeColList[[i]])]
  }
  return(labSamp)
}