Prep_geneInfo <- function(tmp, schemaTmp){
  i <- which(schemaList %in% schemaTmp)
  genesTmp <- subtypeGeneList[[i]]
  idxGene <- which(genesTmp$geneSymbol %in% tmp$featInfo)
  geneInfo <- genesTmp[idxGene,]
  return(geneInfo)
}