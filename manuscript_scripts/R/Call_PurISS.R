Call_PurISS <- function(dataSet, version = "71322") {

  if(version == "71322") {
    load("PurISS/VISIUM_Moff_Dijk_CPTAC_TCGA_Classifer_cut05_alpha05_71322.RData")
    classifier = fitted1$classifier2
  } else if(version == "81822") {
    load("PurISS/VISIUM_Moff_Dijk_CPTAC_TCGA_Classifer_cut05_alpha05_81822.RData")
    classifier = fitted1$classifier2
  } else if(version == "JulyTop10") {
    load("../classifiers/VISIUM_Moff_Dijk_CPTAC_TCGA_Classifer_cut05_alpha05_120522_JulyTop10.RData")
    classifier = fitted1$classifier2
  } else if(version == "Novextreme") {
    load("../classifiers/VISIUM_Moff_Dijk_CPTAC_TCGA_Classifer_cut05_alpha05_120522_Novextreme.RData")
    classifier = fitted1$classifier2
  } else if(version == "Nov") {
    load("../classifiers/VISIUM_Moff_Dijk_CPTAC_TCGA_Classifer_cut05_alpha01_121422_Nov.RData")
    classifier = fitted_alpha01$classifier2
  } else if(version == "012023") {
    load("../classifiers/VISIUM_Moff_Dijk_CPTAC_TCGA_Classifer_cut05_alpha095_012023.Rdata")
    classifier = fitted_alpha095_cut05$classifier2
  } else if(version == "012023_noMoffitt") {
    load("../classifiers/VISIUM_Dijk_CPTAC_TCGA_Classifer_cut25_alpha25_012023.Rdata")
    classifier = fitted_alpha25_cut25_nomoff$classifier2
  } else if(version == "final") {
    load("../classifiers/final_classifier.Rdata")
    classifier = final_classifier$classifier2
  } else {
    stop("version not found!")
  }
  
  source("../R/PurISS/Apply_PurISS.R")
  
  rownames(dataSet$ex) <- make.names(dataSet$featInfo$SYMBOL,unique = TRUE)
  predictions <- Apply_PurISS(data = dataSet$ex, classifier = classifier)
  names(predictions) <- c(paste("PurISS.prob", version, sep = "."),
                          paste("PurISS", version, sep = "."),
                          paste("PurISS_graded", version, sep = ".") )
  
  dataSet$sampInfo <- cbind(dataSet$sampInfo, predictions)
  return(dataSet)
}
