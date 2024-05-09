Call_PurIST <- function(dataSet) {
  rownames(dataSet$ex) <- make.names(dataSet$featInfo$SYMBOL,unique = TRUE)
  classifier = classifs[[1]]
  predictions = apply_classifier(data = dataSet$ex, classifier = classifier)
  
  dataSet$sampInfo$PurIST.prob <- predictions$Pred_prob_basal
  dataSet$sampInfo$PurIST <- predictions$Subtype
  dataSet$sampInfo$PurIST <- gsub("basal-like","Basal-like",dataSet$sampInfo$PurIST)
  dataSet$sampInfo$PurIST <- gsub("classical","Classical",dataSet$sampInfo$PurIST)
  dataSet$sampInfo$PurIST_graded <- predictions$Subtype_graded
  dataSet$sampInfo$PurIST_graded <- gsub("lean basal-like","Lean Basal-like",dataSet$sampInfo$PurIST_graded)
  dataSet$sampInfo$PurIST_graded <- gsub("lean classical","Lean Classical",dataSet$sampInfo$PurIST_graded)
  dataSet$sampInfo$PurIST_graded <- gsub("likely basal-like","Likely Basal-like",dataSet$sampInfo$PurIST_graded)
  dataSet$sampInfo$PurIST_graded <- gsub("likely classical","Likely Classical",dataSet$sampInfo$PurIST_graded)
  dataSet$sampInfo$PurIST_graded <- gsub("strong basal-like","Strong Basal-like",dataSet$sampInfo$PurIST_graded)
  dataSet$sampInfo$PurIST_graded <- gsub("strong classical","Strong Classical",dataSet$sampInfo$PurIST_graded)
  return(dataSet)
}