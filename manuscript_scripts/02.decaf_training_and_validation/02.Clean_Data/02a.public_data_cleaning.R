#####################################################################################

##                             Data Cleaning Script                      ############

####################################################################################




######################################################################################
############## Step 1: Create Functions for Loading Data ########################
######################################################################################




#### Load Raw Subtype Calls ####
Load_data_subtype <- function(rDataName) {
  
 dataSet =  readRDS(sprintf("data/data_DeCAF/%s.caf_subtype.rds", rDataName))
  

 
  ## Change myCAF -> permCAF & iCAF -> restCAF
 
 dataSet$Subtype$SCISSORS_CAF_K2_top25.vK[dataSet$Subtype$SCISSORS_CAF_K2_top25.vK == "myCAF"] = "permCAF"
 dataSet$Subtype$SCISSORS_CAF_K2_top25.vK[dataSet$Subtype$SCISSORS_CAF_K2_top25.vK == "Mixed.myCAF"] = "Mixed.permCAF"
 dataSet$Subtype$SCISSORS_CAF_K2_top25.vK[dataSet$Subtype$SCISSORS_CAF_K2_top25.vK == "Mixed.iCAF"] = "Mixed.restCAF"
  dataSet$Subtype$SCISSORS_CAF_K2_top25.vK[dataSet$Subtype$SCISSORS_CAF_K2_top25.vK == "iCAF"] = "restCAF"
  dataSet$Subtype$CAF = dataSet$Subtype$SCISSORS_CAF_K2_top25.vK
  dataSet$Subtype$SCISSORS_CAF_K2_top25.K2[dataSet$Subtype$SCISSORS_CAF_K2_top25.K2 == "myCAF"] = "permCAF"
  dataSet$Subtype$SCISSORS_CAF_K2_top25.K2[dataSet$Subtype$SCISSORS_CAF_K2_top25.K2 == "iCAF"] = "restCAF"
  
  ## Change Mixed permCAF -> permCAF
  dataSet$Subtype$CAF[dataSet$Subtype$CAF == "Mixed"] = NA
  dataSet$Subtype$CAF[dataSet$Subtype$CAF == "Absent"] = NA
  dataSet$Subtype$CAF[dataSet$Subtype$CAF == "Mixed.permCAF"] = "permCAF"
  dataSet$Subtype$CAF[dataSet$Subtype$CAF == "Mixed.restCAF"] = "restCAF"
  dataSet$Subtype$CAF_K4 = dataSet$Subtype$SCISSORS_CAF_K2_top25.vK
  
  ## remove Absent Call 
  dataSet$Subtype$CAF_ALL = dataSet$Subtype$CAF_K4
  
  dataSet$Subtype$CAF_K4[dataSet$Subtype$CAF_K4 == "Absent"] = NA
  
  dataSet$Subtype$CAF.classifier = dataSet$Subtype$SCISSORS_CAF_K2_top25.vK.classifier
  
  return(dataSet)
}


########################################


#### Load Raw Survival Data ####
Load_data_survival <- function(rDataName) {
  
 dataSet =  readRDS(sprintf("data/data_DeCAF/%s.survival_data.rds", rDataName))
  
  
  return(dataSet)
}




##################################

#### Load Raw Gene and Clinical Data ####
Load_data_gene <- function(rDataName, type = c("array", "seq")) {
  
  
  dat = try(readRDS(sprintf("data/data_DeCAF/%s.rds", rDataName)))
  if(is.character(dat))  {dat = try(readRDS(sprintf("data/data_DeCAF/%s_gencode.rds", rDataName)))}
  if(is.character(dat)) { 
    name_to_split = rDataName
    split_name = strsplit(name_to_split, "\\.")[[1]]
    dat = try(readRDS(
      sprintf("data/data_DeCAF/%s_gencode.%s.rds", split_name[1], split_name[2])))
    
  }
  
  
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
Generate_Survival_Data <- function(survSet, subtypeSet, DeCAF_calls, save, rDataName) {
  
  survDat <- data.frame(ID = survSet$sampID,
                        time = as.numeric(survSet$time),
                        status = as.numeric(survSet$event),
                        surv_whitelist = survSet$whitelist,
                        SCISSORS_CAF = subtypeSet$Subtype$CAF,
                        SCISSORS_CAF_K4 = subtypeSet$Subtype$CAF_K4,
                        SCISSORS_CAF_ALL = subtypeSet$Subtype$CAF_ALL,
                        SCISSORS_CAF_K2 = subtypeSet$Subtype$SCISSORS_CAF_K2_top25.K2,
                        MS_K2 = subtypeSet$Subtype$MS,
                        Elyada_CAF = subtypeSet$Subtype$Elyada_CAF,
                        PurIST = subtypeSet$Subtype$PurIST,
                        Maurer =  subtypeSet$Subtype$Maurer,
                        Puleo = subtypeSet$Subtype$Puleo,
                        DeCAF = DeCAF_calls$Stroma_Subtype,
                        DeCAF_prob = DeCAF_calls$Pred_prob_permCAF,
                        stringsAsFactors = FALSE)
  survDat <- survDat[survSet$whitelist == TRUE,]
  
  if(save == TRUE){
    filename = paste("data/clean_data/survival/", rDataName, "_survival.csv", sep = "")
    write.csv(survDat,filename, row.names = FALSE)
  }
  return(survDat)
}

#################################


## Generate Cleaned Survival data
Generate_SuppTable_Data <- function(survSet, subtypeSet, DeCAF_calls, rDataName) {
  
  SuppTable_Data <- data.frame(Dataset = rDataName,
                               ID = subtypeSet$Subtype$sampID,
                               PDAC_Pi_whitelist = !is.na(subtypeSet$Subtype$CAF_ALL),
                               DeCAF = DeCAF_calls$Stroma_Subtype,
                               DeCAF_graded = DeCAF_calls$Stroma_Subtype_graded,
                               DeCAF_prob = DeCAF_calls$Pred_prob_permCAF,
                               PurIST = subtypeSet$Subtype$PurIST,
                               PurIST_graded = subtypeSet$Subtype$PurIST_graded,
                               PurIST_prob = subtypeSet$Subtype$PurIST.prob,
                               SCISSORS_K2 = subtypeSet$Subtype$SCISSORS_CAF_K2_top25.K2,
                               SCISSORS_vK = subtypeSet$Subtype$SCISSORS_CAF_K2_top25.vK,
                               SCISSORS_combined = subtypeSet$Subtype$CAF,
                               SCISSORS_labels_whitelist = !is.na(subtypeSet$Subtype$CAF) & subtypeSet$Subtype$SCISSORS_CAF_K2_top25.vK.classifier,
                               Elyada = subtypeSet$Subtype$Elyada_CAF,
                               Moffitt_stroma = subtypeSet$Subtype$MS,
                               Maurer =  subtypeSet$Subtype$Maurer,
                               Puleo = subtypeSet$Subtype$Puleo)
  if(class(survSet) != "data.frame"){
    SuppTable_Data$Survival_analysis_whitelist = F
    SuppTable_Data$OS_time = NA
    SuppTable_Data$OS_event = NA
  }
  if(class(survSet) == "data.frame"){
    SuppTable_Data$Survival_analysis_whitelist = survSet$whitelist
    SuppTable_Data$OS_time = as.numeric(survSet$time)
    SuppTable_Data$OS_event = as.numeric(survSet$event)
  }


  return(SuppTable_Data)
}



################################



## Generate Cleaned Clinical Data Function
Generate_Clinical_Data <- function(dataSet, subtypeSet, DeCAF_calls, save, rDataName, survival, survSet) {
  
  ## Get Age
  Age = try(as.numeric(dataSet$sampInfo$Age))
  if(length(Age) == 0){
    Age = try(as.numeric(dataSet$sampInfo$age))
  }
  if(length(Age) == 0){
    Age = try(as.numeric(dataSet$sampInfo$Age.at.initial.pathologic.diagnosis))
  }
  if(length(Age) == 0){
    Age = try( as.numeric(dataSet$sampInfo$Age.at.Diagnosis.in.Years) )
  }
  if(length(Age) == 0){
    Age = try( as.numeric(dataSet$sampInfo$age_at_diagnosis) )
  }
  if(length(Age) == 0){
    Age = try( as.numeric(dataSet$sampInfo$Age_at_diagnosis) )
  }
  if(length(Age) == 0){
    Age = try( as.numeric(dataSet$sampInfo$`age (year)`))
  } 
  if(length(Age) == 0){
    Age = NA
  }
  
  
  ## Get Sex
  Sex = try(as.character(dataSet$sampInfo$Gender))
  if(length(Sex) == 0){
    Sex = try(as.character(dataSet$sampInfo$Sex))
  }
  if(length(Sex) == 0){
    Sex = try(as.character(dataSet$sampInfo$sex ))
  }
  if(length(Sex) == 0){
    Sex = try(as.character(dataSet$sampInfo$gender ))
  }
  if(length(Sex) == 0){
    Sex = NA
  }
  
  
  ## Get Race
  Race = try(as.character(dataSet$sampInfo$race))
  if(length(Race) == 0){
    Race = try(as.character(dataSet$sampInfo$Race))
  }
  if(length(Race) == 0){
    Race = try( as.character(dataSet$sampInfo$Ethnicity ))
  }
  if(length(Race) == 0){
    Race = try( as.character(dataSet$sampInfo$ethnicity ) )
  }
  if(length(Race) == 0){
    Race = NA
  }
  
  ## Get Margin
  Margin = try(as.character(dataSet$sampInfo$Margin))
  if(length(Margin) == 0){
    Margin = try(as.character(dataSet$sampInfo$margin))
  }
  if(length(Margin) == 0){
    Margin = try( as.character(dataSet$sampInfo$Residual.tumor ))
  }
  if(length(Margin) == 0){
    Margin = try( as.character(dataSet$sampInfo$Characteristics.resection.margin. ) )
  } 
  if(length(Margin) == 0){
    Margin = NA
  } 
  
  
  ## Generate Data Frame
  clinicalDat <- data.frame(Sex =  Sex,
                            Race = Race,
                            Age = Age,
                            Margin = Margin,
                            SCISSORS_CAF = subtypeSet$Subtype$CAF,
                            SCISSORS_CAF_K4 = subtypeSet$Subtype$CAF_K4,
                            SCISSORS_CAF_ALL = subtypeSet$Subtype$CAF_ALL,
                            MS_K2 = subtypeSet$Subtype$MS,
                            Elyada_CAF = subtypeSet$Subtype$Elyada_CAF,
                            PurIST = subtypeSet$Subtype$PurIST,
                            Maurer =  subtypeSet$Subtype$Maurer,
                            Puleo = subtypeSet$Subtype$Puleo,
                            CAF_whitelist =  subtypeSet$Subtype$CAF.classifier,
                            stringsAsFactors = FALSE)
  
  if(survival == TRUE){
    clinicalDat$surv_whitelist = survSet$whitelist
  } else {
    clinicalDat$surv_whitelist = NA
  }
  
  ## Standardize Categories
  clinicalDat <- clinicalDat[clinicalDat$CAF_whitelist == TRUE,]
  
  clinicalDat$Sex <- gsub("female","F",clinicalDat$Sex)
  clinicalDat$Sex <- gsub("male","M",clinicalDat$Sex)
  clinicalDat$Sex <- gsub("FeM","F",clinicalDat$Sex)
  clinicalDat$Sex <- gsub("Male","M",clinicalDat$Sex)
  clinicalDat$Sex <- gsub("NaN",NA,clinicalDat$Sex)
  clinicalDat$Race <- gsub("white","White",clinicalDat$Race)
  clinicalDat$Race <- gsub("White/Caucasian","White",clinicalDat$Race)
  clinicalDat$Race <- gsub("Black/African","AA",clinicalDat$Race)
  clinicalDat$Race <- gsub("Asian, Caucasian","Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("Not documented","Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("black or african american","AA",clinicalDat$Race)
  clinicalDat$Race <- gsub("Black or African American","AA",clinicalDat$Race)
  clinicalDat$Race <- gsub("Unknown or Not Reported",NA,clinicalDat$Race)
  clinicalDat$Race <- gsub("Unknown",NA,clinicalDat$Race)
  clinicalDat$Race <- gsub("African American","AA",clinicalDat$Race)
  clinicalDat$Race <- gsub("Caucasian/CaucOther","Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("Caucasian_Non-Hispanic","White",clinicalDat$Race)
  clinicalDat$Race <- gsub("Caucasian_Latino","White",clinicalDat$Race)
  clinicalDat$Race <- gsub("CaucOther","Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("Hispanic_Latino","Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("Hispanic","Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("asian","Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("Asian","Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("Cauc","White",clinicalDat$Race)
  clinicalDat$Race <- gsub("WhiteOther", "Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("Other_Latino", "Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("Other, White", "Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("Other, White", "Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("WhiteOther", "Other",clinicalDat$Race)
  
  
  
  clinicalDat$Race <- gsub("Black","AA",clinicalDat$Race)
  clinicalDat$Race <- gsub("CaucOther","Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("Hispanic","Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("Hispanic_Latino","Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("Other, Caucasian/CaucOther","Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("White_Latino","White",clinicalDat$Race)
  clinicalDat$Race <- gsub("White_Non-Other","Other",clinicalDat$Race)
  clinicalDat$Race <- gsub("Caucasian", "White",clinicalDat$Race)
  clinicalDat$Race <- gsub("caucasian", "White",clinicalDat$Race)
  
  
  
  clinicalDat$Margin <- gsub("negative","R0",clinicalDat$Margin)
  clinicalDat$Margin <- gsub("close","R0",clinicalDat$Margin)
  clinicalDat$Margin <- gsub("positive","R1",clinicalDat$Margin)
  clinicalDat$Margin <- gsub("r0","R0",clinicalDat$Margin)
  clinicalDat$Margin <- gsub("r1","R1",clinicalDat$Margin)
  clinicalDat$Margin <- gsub("resection margin ","",clinicalDat$Margin)
  clinicalDat$Margin <- gsub("not available",NA,clinicalDat$Margin)
  
  clinicalDat$Age <- round(clinicalDat$Age)
  
  
  if(save == TRUE){
    filename = paste("data/clean_data/clinical/", rDataName, "_clinical.csv", sep = "")
    write.csv(clinicalDat,filename, row.names = FALSE)
  }
  return(clinicalDat)
}



## Generate Cleaned Subtype data
Generate_Subtype_Data <- function(subtypeSet, DeCAF_calls, save, rDataName) {
  
  subtype <- data.frame(ID = subtypeSet$Subtype$sampID,
                        SCISSORS_CAF = subtypeSet$Subtype$CAF,
                        SCISSORS_CAF_K4 = subtypeSet$Subtype$CAF_K4,
                        SCISSORS_CAF_ALL = subtypeSet$Subtype$CAF_ALL,
                        SCISSORS_CAF_K2 = subtypeSet$Subtype$SCISSORS_CAF_K2_top25.K2,
                        MS_K2 = subtypeSet$Subtype$MS,
                        Elyada_CAF = subtypeSet$Subtype$Elyada_CAF,
                        PurIST = subtypeSet$Subtype$PurIST,
                        Maurer =  subtypeSet$Subtype$Maurer,
                        Puleo = subtypeSet$Subtype$Puleo,
                        DeCAF = DeCAF_calls$Stroma_Subtype,
                        CAF_whitelist = subtypeSet$Subtype$CAF.classifier,
                        stringsAsFactors = FALSE)
  subtype <- subtype[subtype$CAF_whitelist == TRUE,]
  
  if(save == TRUE){
    filename = paste("data/clean_data/subtype/", rDataName, "_subtype.csv", sep = "")
    write.csv(subtype,filename, row.names = FALSE)
  }
  return(subtype)
}


################################



## Extract Classifier Genes
Extract_Classifier_Genes <- function(data, classifier, subtypeSet, save, rDataName) {
 
  rows = which(rownames(data) %in% classifier$TSPs[classifier$fit$beta[-1] != 0,])
  white_list = subtypeSet$Subtype$CAF.classifier
  if(sum(is.na(white_list)) == nrow(subtypeSet$Subtype)){
    white_list = rep(T, nrow(subtypeSet$Subtype))
  }
  data = data[rows,white_list]
  tsp = match(as.vector(t(classifier$TSPs[classifier$fit$beta[-1]!=0,])),rownames(data))
  data = data[tsp,]
  
  
  if(save == TRUE){
    filename = paste("data/clean_data/classifier_gene_exp/", rDataName, "_classifier_gene_expression.csv", sep = "")
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
training = unique(final_DeCAF_classifier$anno$study)




######################################################################################
######## Step 3: Clean Survival, Subtype, and Clinical Raw Data #######################
######################################################################################




## Names of Data
DataNames <- c("CPTAC", "Dijk", "Grunwald", "Hayashi",
               "Linehan", "Moffitt_GEO_array", "Olive", "PACA_AU_array", "PACA_AU_seq",
               "Puleo_array", "TCGA_PAAD")
## Data Types
data_type = c("seq", "seq", "seq", "seq",
              "seq", "array", "seq", "array", "seq",
              "array", "seq")
## Is there survival info?
survival_info = c(T, T, T, F,
              T, T, F, T, T,
              T, T)



## Loop Through All Data
for(i in 1:length(DataNames)){
  
  print(paste(DataNames[i], data_type[i], "survival info", survival_info[i], sep = " - "))
  
  ## Load Data
  rDataName <- DataNames[i]
  subtypeSet <- Load_data_subtype(rDataName)
  dataSet <- Load_data_gene(rDataName, type = data_type[i] )
  print("loading done")
  
  ## Generate DeCAF Calls
  DeCAF_calls <- evaluate_classifier(dataSet$expression, classifier)
  
  ## Clean and Save Survival Data if applicable
  if(survival_info[i]){
    survSet <- Load_data_survival(rDataName)
    survDat <- Generate_Survival_Data(survSet, subtypeSet, DeCAF_calls, save = TRUE, rDataName)
    survDat$dataSet <- rDataName
    if(i == 1){
      survDatCmb <- survDat
    } else {
      survDatCmb <- rbind(survDatCmb,survDat)
    }
    
  }
  if(!survival_info[i]){
    survSet = NA
  }
  
  ## Save for Supplemental Table
  SuppTable_Data_study = Generate_SuppTable_Data(survSet, subtypeSet, DeCAF_calls, rDataName)
  
  ## Clean and Save Data
  clinicalDat <- Generate_Clinical_Data(dataSet, subtypeSet, DeCAF_calls, save = TRUE, 
                                        rDataName, survival = survival_info[i], survSet)
  subtypeDat <- Generate_Subtype_Data(subtypeSet, DeCAF_calls, save = TRUE, rDataName)
  classif_exp <- Extract_Classifier_Genes(dataSet$expression, classifier, subtypeSet, save = T, rDataName)
  
  
  ## Add Study Name
  clinicalDat$dataSet <- rDataName
  subtypeDat$dataSet <- rDataName
  DeCAF_calls$dataSet <- rDataName
  
  ## Only Subtype True DeCAF Calls
  DeCAF_calls = DeCAF_calls[subtypeSet$Subtype$CAF.classifier,]
  
  ## Combine With Other Studies
  if(i == 1){
    clinicalCmb <- clinicalDat
    subtypeDatCmb <- subtypeDat
    classif_expCmb <- classif_exp
    DeCAF_callsCmb <- DeCAF_calls
    SuppTable_Data <- SuppTable_Data_study
  } else {
    clinicalCmb <- rbind(clinicalCmb,clinicalDat)
    subtypeDatCmb <- rbind(subtypeDatCmb, subtypeDat)
    classif_expCmb <- cbind(classif_expCmb, classif_exp)
    DeCAF_callsCmb <- rbind(DeCAF_callsCmb, DeCAF_calls)
    SuppTable_Data <- rbind(SuppTable_Data, SuppTable_Data_study)
  }
  
  
  
  
}




## Save Combined Data
write.csv(survDatCmb,"data/clean_data/combined_data/combined_survival_DeCAF.csv", row.names = FALSE)
write.csv(subtypeDatCmb,"data/clean_data/combined_data/combined_subtype_DeCAF.csv", row.names = FALSE)
write.csv(classif_expCmb,"data/clean_data/combined_data/combined_classifier_gene_exp_DeCAF.csv")
write.csv(DeCAF_callsCmb,"data/clean_data/combined_data/combined_classifier_calls_DeCAF.csv")

## Change Survival of Duplicated ID's to Individual Level Only
idDup <- SuppTable_Data$ID[(duplicated(SuppTable_Data$ID))]
idDup_Paca_array <- which((SuppTable_Data$ID %in% idDup) & (SuppTable_Data$dataSet == "PACA_AU_array"))
SuppTable_Data$Survival_analysis_whitelist[idDup_Paca_array] <- "Individual Study Only - Not Combined"
write.xlsx(SuppTable_Data, "data/Supp_Table/Supplementary Table S5 - clinical and molecular data - public.xlsx", row.names = F)



