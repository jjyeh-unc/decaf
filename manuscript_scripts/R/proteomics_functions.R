convertToMouseSymbol <- function(symbol){
  original <- symbol
  converted <- paste(toupper(substr(original, 1, 1)), 
                     substr(tolower(original), 2, nchar(tolower(original))), sep="")
  return(converted)
}

loadKinaseDictionary <- function(){
  kinaseDictionary <- read.table(file = "./Proteomics_utitilies/Updated_Kinase_List_For_150929_Uniprot_FASTA.csv",
                                 header = TRUE,
                                 sep = ",",
                                 stringsAsFactors = FALSE)
  
    mouse <- cbind(kinaseDictionary[,c(1,2,4,5)])
    mouse$Gene.Symbol <- paste(toupper(substr(mouse$Gene.Symbol, 1, 1)), 
                               substr(tolower(mouse$Gene.Symbol), 2, nchar(tolower(mouse$Gene.Symbol))), sep="")
    mouse$Uniprot.Kinase.Name <- paste(toupper(substr(mouse$Uniprot.Kinase.Name, 1, 1)), 
                                       substr(tolower(mouse$Uniprot.Kinase.Name), 2, nchar(tolower(mouse$Uniprot.Kinase.Name))), sep="")
    names(mouse)[3] <- "Protein.IDs"
    mouse$species <- "mouse"
    
    human <- cbind(kinaseDictionary[,c(1,2,3,5)])
    names(human)[3] <- "Protein.IDs"
    human$species <- "human"
    
    
    kinaseDictionary <- rbind(human,mouse)
    kinaseDictionary <- subset(kinaseDictionary,!(kinaseDictionary$Protein.IDs==""))
    names(kinaseDictionary)[2] <- "Kinase.Name"
          
    return(kinaseDictionary)
}

getCleanKinase <- function(queryDF, queryCol = "Gene.names", queryType = "geneSymbol", geneAsRownames = TRUE){
  query_df <- queryDF
  query_col <- queryCol
  query_type <- queryType
  
  kinaseDictionary <- loadKinaseDictionary()
  if(query_type == "geneSymbol"){
    kinaseList <- unlist(kinaseDictionary$Gene.Symbol)
  }  else if(query_type == "proteinID"){
    kinaseList <- unlist(kinaseDictionary$Protein.IDs)
  }  else if(query_type == "proteinSymbol"){
    kinaseList <- unlist(kinaseDictionary$Kinase.Name)
  }
  
  checkKinaseHit <- function(queryList){
    queryList <- strsplit(as.character(queryList), split=';', fixed=TRUE)[[1]]
    hit_cnt <- as.numeric(sum(queryList %in% kinaseList))
    if(hit_cnt == 1){
      matched_kinase <- as.matrix(kinaseDictionary$Kinase.Name[kinaseList %in% queryList])
      return(matched_kinase)
    } else {
      return("DROP")
    }
  }
        
  Kinase.names <- as.matrix(as.character(sapply(query_df[ ,query_col], function(x) checkKinaseHit(x)) ))
  
  res_df <- query_df[,!(names(query_df) %in% query_col)]
  res_df <- cbind(Kinase.names, res_df)
  names(res_df)[1] <- "Kinase.names"
  res_df <- res_df[!(res_df$Kinase.names == "DROP"),]
  
  if(geneAsRownames) {
    row.names(res_df) <- make.names(kinaseDictionary$Gene.Symbol[match(res_df$Kinase.names, kinaseDictionary$Kinase.Name)], unique=TRUE)
  }
  
  # sort the kinases
  res_df <- res_df[order(res_df$Kinase.name), ]
  species <- res_df$Kinase.name == toupper(res_df$Kinase.name)
  res_df <- res_df[order(-species), ]
  return(res_df)
}


readMaxQuantProtein <- function(filename = "./maxquant/proteinGroups.txt", 
                                LFQ = TRUE,
                                dropREV = TRUE,
                                dropCON = TRUE,
                                noMissing = TRUE,
                                proteinType = "kinaseHit"|"proteinHit"|"allProteins"){
  rawdata <- read.table(file = filename,
                        header = TRUE,
                        sep = "\t")

  output <- rawdata[ ,c("Majority.protein.IDs", "Gene.names"), drop=FALSE]
  names(output)[1] <- "Protein.IDs"
  
  if(LFQ){
    intensity_cols <- grep("LFQ.intensity",names(rawdata))
  } else {
    intensity_tmp <- names(rawdata)[grep("Reporter.intensity", names(rawdata))]
    ele_num <- sapply(intensity_tmp, function(x) length(strsplit(x, split='.', fixed=TRUE)[[1]]))
    intensity_cols <- which(names(rawdata) %in% intensity_tmp[which(ele_num==3)])
  }
  
  output <- cbind(output,rawdata[,intensity_cols,drop=FALSE])
  
  if(dropREV){
    reversed <- grep("REV__",output$Protein.IDs)
    output <- output[-reversed,]
  }
  
  if(dropCON){
    CON <- grep("CON__",output$Protein.IDs)
    output <- output[-CON,]
  }
  
  if(noMissing){
    expr <- output[ , -which(names(output) %in% c("Protein.IDs", "Gene.names") )]
    feature <- which(apply(expr,1,FUN = function(x) sum(x==0))<=0)
    expr <- log10(expr[feature, ]+1)
    
    output <- output[feature, ]
    output <- cbind(output[ ,c("Protein.IDs","Gene.names")],expr)
  }
  
  if(proteinType == "kinaseHit"){
    output <- output[ ,!(names(output) %in% "Gene.names")]
    output <- getCleanKinase(queryDF = output, queryCol = "Protein.IDs", queryType = "proteinID", geneAsRownames = TRUE)
  } else if(proteinType == "proteinHit"){
    ambiguous <- grep(";",output$Protein.IDs)
    output <- output[-ambiguous,]
    output <- output[ ,!(names(output) %in% "Protein.IDs")]
    rownames(output) <- make.names(output[, 1], unique=TRUE)
  } 
  
  return(output)
}

readProteinPilotProtein <- function(filename = "./ProteinPilot/ProteinPilot_Yeh15_human_ProteinSummary.txt",
                                    dropREV = TRUE,
                                    proteinType = "kinaseHit"|"proteinHit"|"allProteins"){
  rawdata <- read.table(file = filename,
                        header = TRUE,
                        sep = "\t",
                        quote = "",
                        check.names=FALSE)
  
  rawdata[sapply(rawdata, function(x) all(is.na(x)))] <- NULL
  
  if(dropREV){
    reversed <- grep(x = as.character(rawdata[,c("Accession")]),
                     fixed = TRUE,
                     pattern = "RRRRRsp")
    rawdata <- rawdata[-reversed,]
  }
  
  ratio_cols <- names(rawdata)[grep(":",names(rawdata))] 
  expr <- rawdata[ ,(names(rawdata) %in% ratio_cols[-grep(" ", ratio_cols)])]
  features <- which(apply(expr,1,FUN = function(x) sum(x==0))<=0)
  
  output <- cbind(rawdata[ ,c("Accession","Species")], expr)[features, ]
  proteinInfo <- strsplit(as.character(output[,c("Accession")]),
                          fixed=TRUE, split="|")
  proteinInfo <- strsplit(sapply(proteinInfo, function(x) x[3]),
                          fixed=TRUE,split="_")
  Protein.symbol <- as.character(sapply(proteinInfo, function(x) x[1]))
  output <- cbind(Protein.symbol, output)
  output <- output[, !(names(output) %in% c("Accession"))]
  output$Protein.symbol <- as.character(output$Protein.symbol)
  output$Protein.symbol[which(output$Species %in% "MOUSE")] <- as.character(convertToMouseSymbol(output$Protein.symbol[output$Species == "MOUSE"]))
 
  output <- output[ , -which(names(output) %in% c("Species") )]  
  
  if(proteinType == "kinaseHit"){
    output <- getCleanKinase(queryDF = output, queryCol = "Protein.symbol", queryType = "proteinSymbol", geneAsRownames = TRUE)
  } else if(proteinType == "proteinHit"){
   # ambiguous <- grep(";",output$Protein.symbol)
    #output <- output[-ambiguous,]
    rownames(output) <- make.names(output[, 1], unique=TRUE)
  } 
  
  return(output)
}

readProteinPilotPeptide <- function(filename = "./ProteinPilot/ProteinPilot_Yeh15_human_PeptideSummary.txt",
                                    denominator = "114",
                                    dropREV = TRUE,
                                    dropAmbiguous = TRUE,
                                    proteinType = "kinaseHit"|"proteinHit"){
  rawdata <- read.table(file = filename,
                        header = TRUE,
                        sep = "\t",
                        quote = "",
                        check.names=FALSE)
  
  rawdata[sapply(rawdata, function(x) all(is.na(x)))] <- NULL
  
  if(dropREV){
    reversed <- grep(x = as.character(rawdata[,c("Accessions")]),
                     fixed = TRUE,
                     pattern = "RRRRRsp")
    rawdata <- rawdata[-reversed,]
  }
  
  output <- cbind(rawdata[ ,c("Accessions","Used")])

  # read in peptide level file and summarize to protein ratio using perl
  expr_area <- rawdata[grep("Area",names(rawdata))] 
  expr_area[is.na(expr_area)] <- 0
  name_rf <- paste("Area ",denominator,sep="")
  expr_rf <- expr_area[ names(expr_area) %in% name_rf]
  expr_tx <- expr_area[ !(names(expr_area) %in% name_rf)]
  
  for(i in names(expr_tx)) {
    ratio_name <- paste(sub("Area ", "", i),":",denominator,sep="")
    
    tmp_tx <- expr_tx[i] 
    expr_tmp <- cbind(expr_rf, tmp_tx)
    good_rows <- which(apply(expr_tmp,1,FUN = function(x) sum(x==0))<=0)
    
    ratio_tmp <- as.numeric(expr_tmp[good_rows, 2]/expr_tmp[good_rows, 1])
    ratio_norm <- ratio_tmp/median(ratio_tmp)
    output[good_rows ,paste(ratio_name)] <- as.numeric(ratio_norm)
    output[!good_rows,paste(ratio_name)] <- NA
    
    rm(ratio_name, tmp_tx, expr_tmp, ratio_tmp, ratio_norm)
  }
  
  if(dropAmbiguous){
    ambiguous <- grep(x = as.character(output[,c("Accessions")]),
                     fixed = TRUE,
                     pattern = ";")
    output <- output[-ambiguous,]
  }
  
  peptideInfo <- strsplit(as.character(output[,c("Accessions")]),
                          fixed=TRUE, split="|")
  peptideInfo <- strsplit(sapply(peptideInfo, function(x) x[3]),
                          fixed=TRUE,split="_")
  Protein.symbol <- as.character(sapply(peptideInfo, function(x) x[1]))
  Species <- as.character(sapply(peptideInfo, function(x) x[2]))
  output <- cbind(Protein.symbol, Species, output)
  output <- output[, !(names(output) %in% c("Accessions"))]
  output$Protein.symbol <- as.character(output$Protein.symbol)
  output$Protein.symbol[which(output$Species %in% "MOUSE")] <- as.character(convertToMouseSymbol(output$Protein.symbol[output$Species == "MOUSE"]))
  
  output <- output[ , -which(names(output) %in% c("Species") )]  
  
  output <- output[output$Used != 0, ] 
  output <- output[ , -which(names(output) %in% c("Used") )]  
  
  output <- aggregate(.~Protein.symbol, output, median, na.action=na.omit)
  
  if(proteinType == "kinaseHit"){
    output <- getCleanKinase(queryDF = output, queryCol = "Protein.symbol", queryType = "proteinSymbol", geneAsRownames = TRUE)
  } else if(proteinType == "proteinHit"){
    rownames(output) <- make.names(output[, 1], unique=TRUE)
  } 
  
  
  return(output)
}
