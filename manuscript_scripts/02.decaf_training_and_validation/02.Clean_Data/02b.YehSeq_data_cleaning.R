#####################################################################################

##                             Data Cleaning Script                      ############

####################################################################################




######################################################################################
############## Step 1: Create Functions for Loading Data ########################
######################################################################################



#### Load Raw Survival Data ####
Load_yehseqdata_survival <- function() {
  
  dataSet =  readRDS("data/YehSeq_Pi/yehseq_pdac_pi.salmon_refseq.no_tpl.rds")
  
  
  return(dataSet)
}




##################################

#### Load Raw Gene and Clinical Data ####
Load_yehseqdata_gene <- function() {
  
  type = "seq"
  dat = readRDS("data/YehSeq_Pi/yehseq_pdac_pi.salmon_refseq.no_tpl.rds")
  
  
  # process expression information
  dat2 = dat$ex
  dat2 = dat2[!is.na(dat$featInfo$SYMBOL),]
  genes = dat$featInfo$SYMBOL[!is.na(dat$featInfo$SYMBOL)]
  if(any(table(genes) > 1)){
    # need to collapse expression
    if(type == "array"){
      dat2 = aggregate(. ~ genes, dat = dat2, FUN = mean)
      rownames(dat2) = dat2[,1]
      dat2 = dat2[,-1]
    }else{
      dat2 = aggregate(. ~ genes, dat = data.frame(dat2), FUN = sum)
      rownames(dat2) = dat2[,1]
      dat2 = dat2[,-1]
      #dat2 = dat2[rowMeans(dat2 > 5) > .2 ,]
    }
  }else{
    rownames(dat2) = dat$featInfo$SYMBOL
  }
  
  colnames(dat2) = colnames(dat$ex)     
  
  ## Resave back to dat
  dat$expression = dat2
  
  
  return(dat)
}



#################################


## Generate Cleaned Survival data
Generate_Survival_Data <- function(dataSet, DeCAF_calls, save) {
  
  survDat <- data.frame(ID = colnames(dataSet$ex),
                        time = dataSet$clinicalDat$`Overall.Survival.(OR.until.death)`,
                        event = dataSet$clinicalDat$Last.Event,
                        neoadj.tx = dataSet$clinicalDat$Neoadj.Tx,
                        neoadj.tx.clean = dataSet$clinicalDat$neo_reg_clean,
                        M = dataSet$clinicalDat$M,
                        surv_whitelist = TRUE,
                        PurIST = dataSet$sampInfo$PurIST,
                        DeCAF = DeCAF_calls$Stroma_Subtype,
                        DeCAF_prob = DeCAF_calls$Pred_prob_permCAF,
                        stringsAsFactors = FALSE)
  
  survDat$event[survDat$event == 500 | 
                  survDat$event ==900 |
                  survDat$event ==600 |
                  survDat$event ==650 ] <- 0
  survDat$event[survDat$event == 2 | 
                  survDat$event ==3] <- 1
  survDat$event = as.numeric(survDat$event)
  
  colnames(survDat)[colnames(survDat) == "event"] = "status"
  survDat$time <- as.numeric(survDat$time)
  
                      
  if(save == TRUE){
    filename = paste("data/clean_data/YehSeq/YehSeq_PDAC_Pi.salmon_refseq.no_tpl.survival.csv", sep = "")
    write.csv(survDat,filename, row.names = FALSE)
  }
  return(survDat)
}
################################





## Extract Classifier Genes
Extract_Classifier_Genes <- function(data, classifier, save) {
  
  rows = which(rownames(data) %in% classifier$TSPs[classifier$fit$beta[-1] != 0,])
  
  data = data[rows,]
  tsp = match(as.vector(t(classifier$TSPs[classifier$fit$beta[-1]!=0,])),rownames(data))
  data = data[tsp,]
  
  
  if(save == TRUE){
    filename = paste("data/clean_data/YehSeq/YehSeq_PDAC_Pi.salmon_refseq.no_tpl.classifier_gene_expression.csv", sep = "")
    write.csv(data,filename, row.names = FALSE)
  }
  return(data)
}



######################################################################################
################ Step 2: Load Classifier #############################################
######################################################################################

## Set working directory
## setwd(PATH TO DeCAF_Classifier_Clean folder)

## Load Classifier
load("DeCAF_classifier/final_DeCAF_classifier.Rdata")
classifier = final_DeCAF_classifier$classifier2

## Load Evaluate Classifier Function 
source("scripts/00.Misc_Functions/evaluate_DeCAF_classifier_function.R")




######################################################################################
######## Step 3: Clean Survival, Subtype, and Clinical Raw Data #######################
######################################################################################



################ Load Salmon Refseq Data


  ## Load YehSeq Gene Data

  dataSet <- Load_yehseqdata_gene()
  print("loading done")
  
  ## Generate DeCAF Calls
  DeCAF_calls <- evaluate_classifier(dataSet$expression, classifier)
  
  survSet <- Load_yehseqdata_survival()
  survDat <- Generate_Survival_Data(survSet, DeCAF_calls, save = TRUE)
  survDat$dataSet <- "YehSeq_PDAC_Pi"
  
  classif_exp <- Extract_Classifier_Genes(dataSet$expression, classifier, save = T)
  DeCAF_calls$dataSet <- "YehSeq_PDAC_Pi"
  
 
  
  write.csv(DeCAF_calls,"data/clean_data/YehSeq/YehSeq_PDAC_Pi.salmon_refseq.no_tpl.classifier_calls.csv")




