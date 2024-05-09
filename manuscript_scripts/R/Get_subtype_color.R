Get_subtype_color <- function(callsTmp, schemaTmp) {
  i <- which(schemaList %in% schemaTmp)
  colTmp <- subtypeColList[[i]]
  if(class(callsTmp) == "numeric") {
    colSamp <- colTmp[callsTmp]
  } else {
    colSamp <- colTmp[match(callsTmp,subtypeList[[i]])]
  }
  return(colSamp)
}