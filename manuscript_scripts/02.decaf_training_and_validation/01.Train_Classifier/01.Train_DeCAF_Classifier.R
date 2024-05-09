
## Load libraries
library(xlsx)
library(survival)
library(survminer)
library(stringr)

####################################################################
####################### LOAD DATA ##################################
####################################################################

# get study names with caf subtypes
study_subtype_file = list.files(path = "data/training_data/",pattern = "*caf_subtype.rds")
study_subtype_names = gsub("\\.caf_subtype.rds","",study_subtype_file)


# for each data set, load the file, corresponding subtypes, and survival data (if applicable)
## Important to load testing data in order to ensure that classifier genes are those present in all available studies
## Makes it more likely that these genes will be present in future studies
for(i in 1:length(study_subtype_names)){
    
    # print study names
    print(study_subtype_names[i])
    
    # read in subtype file
    sub = readRDS(sprintf("data/training_data/%s.caf_subtype.rds", study_subtype_names[i] )) 
    
    
    ## Create DF 
      subDat <- data.frame(
        study = rep(study_subtype_names[i], nrow(sub$Subtype)),
        CAF_K4 = sub$Subtype$SCISSORS_CAF_K2_top25.vK,
        CAF_K2 = sub$Subtype$SCISSORS_CAF_K2_top25.K2,
        Moffitt = sub$Subtype$PurIST,
        Elyada_CAF = sub$Subtype$Elyada_CAF,
        MS_K2 = sub$Subtype$MS,
        sub_whitelist = sub$Subtype$SCISSORS_CAF_K2_top25.vK.classifier,
        sub_whitelist2 = sub$Subtype$SCISSORS_CAF_K2_top25.K2.classifier,
        ID = sub$Subtype$sampID,
        Maurer = sub$Subtype$Maurer,
        Puleo = sub$Subtype$Puleo,
        stringsAsFactors = FALSE
      )
    
    
    # read in expression data file
    dat = try(readRDS(sprintf("data/training_data/%s.rds", study_subtype_names[i])))
    
    
    # process expression information
    dat2 = dat$ex
    dat2 = dat2[!is.na(dat$featInfo$SYMBOL),]
    genes = dat$featInfo$SYMBOL[!is.na(dat$featInfo$SYMBOL)]
    if(any(table(genes) > 1)){
      # need to collapse expression
      if(length(grep("array", study_subtype_names[i]))>0){
        dat2 = aggregate(. ~ genes, dat = dat2, FUN = mean)
        rownames(dat2) = dat2[,1]
        dat2 = dat2[,-1]
      }else{
        dat2 = aggregate(. ~ genes, dat = data.frame(dat2), FUN = sum)
        rownames(dat2) = dat2[,1]
        dat2 = dat2[,-1]
       }
    }else{
      rownames(dat2) = dat$featInfo$SYMBOL
    }
    colnames(dat2) = colnames(dat$ex)     
    
    # save info
    if(i == 1){
      subDatf = subDat
      datf = dat2
    }else{
      subDatf = rbind(subDatf, subDat)
      datf = merge(datf, dat2,by="row.names")
      rownames(datf) = datf[,1]
      datf = datf[,-1]
    }
    
  }
 


####################################################################
################## GENERATE TRAINING FUNCTION ######################
####################################################################

# source helper functions
source("scripts/00.Misc_Functions/classifier_training_functions.R")

## Create classifier function
create.classifier = function(datf, anno, genes0, keep, training0, cut = 0.005, alpha = 0.5){
  keepind = which(anno[,1] %in% keep)
  group = as.character(anno[keepind,1])
  datf = datf[,keepind]
  anno = anno[keepind,]
  
  tumor = which(anno[,2] %in% c("myCAF","iCAF"))
  print(length(tumor))
  
  datf1 = datf[rownames(datf) %in% genes0,tumor]
  anno1 = anno[tumor,]
  group1 = group[tumor]
  s = rep(0, nrow(anno1))
  s[anno1[,2] == "myCAF"] = 1
  s = factor(s)
  
  cnames = colnames(datf1)
  rnames = rownames(datf1)
  datf1 = matrix(as.numeric(as.matrix(datf1)), ncol = ncol(datf1), nrow = nrow(datf1))
  rownames(datf1) = rnames
  colnames(datf1) = cnames
  
  # apply TSP
  library(switchBox)
  datf2 = apply(datf1, 2, rank)
  p = dir = matrix(0, nrow(datf1), length(unique(anno1[,1])))
  for(i in 1:nrow(datf1)){
    for(j in 1:length(unique(anno1[,1]))){
      ok = which(anno1[,1] == unique(anno1[,1])[j])
      p[i,j] = wilcox.test(datf1[i,ok] ~ s[ok])$p.value
      dir[i,j] = mean(datf2[i,ok][s[ok] == 1]) - mean(datf2[i,ok][s[ok] == 0])
    }
  }
  mr = rowMeans((apply(p, 2, rank)))
  top = mr < quantile(mr, cut)
  stop = T
  condition = top & (rowSums(dir < 0) == length(unique(anno1[,1])) | rowSums(dir < 0) == 0)
  
  classifier2 <- SWAP.KTSP.Train(inputMat = datf1[condition,], phenoGroup = s,featureNo=100)
  classifier2 = fitfunc(train_sub =datf1 , classifier = classifier2, class = s, skip = T, alpha = alpha)
  res = create.classif(dat=datf1, classifier=classifier2,dec = dec, labels = s, fit = classifier2$fit)
  trainingPrediction2 <- SWAP.KTSP.Classify(datf1, classifier2)
  trainingPrediction = res$predprob
  print(table(s, trainingPrediction>.5)/length(s))
  return(list(datf1 = datf1, s=s, classifier2 = classifier2, trainingPrediction = trainingPrediction, anno = anno1, group = group1))
}


####################################################################
###################### SET UP DATA FOR TRAINING ####################
####################################################################

#Set up anno - class labels for training
anno = data.frame(study = subDatf$study, subtype = subDatf$CAF_K4)
anno$subtype[anno$subtype == "Mixed.iCAF"] = "iCAF"
anno$subtype[anno$subtype == "Mixed.myCAF"] = "myCAF"
anno$subtype[anno$subtype == "Absent"] = NA
anno$subtype[anno$subtype == "Mixed"] = NA

## Remove white listed samples
anno2 = anno[subDatf$sub_whitelist == T,]  

# setup datf - Gene expression for training and remove white listed samples
datf2 = datf[,subDatf$sub_whitelist == T ]

# setup training sets to keep
sort(table(anno2$study))
training0 = c("Moffitt_GEO_array","Dijk","CPTAC", "TCGA_PAAD")

## setup genes0 - list of genes associated with SCISSORS -> used for classifier
ge = read.xlsx("data/SCISSORS_CAF_marked.xlsx", sheetIndex = 1)
ge = ge[ge$cluster %in% c("myCAF","iCAF"),]
ge = ge[ge$Visium_genes == TRUE,]
genes0 = as.character(ge$gene)

## Remove genes with orf
genes0 = genes0[!grepl("orf", genes0)]


####################################################################
######################### CREATE AND SAVE CLASSIFIER ###############
####################################################################

# Create classifier
fitted_alpha25_cut05_K2 = create.classifier(datf = datf2, anno = anno2, genes0 = genes0, keep = training0, training0 = training0, cut = 0.25, alpha = 0.05)

## Rename as final DeCAF classifer and change labels to be permCAF and restCAF
final_DeCAF_classifier = fitted_alpha25_cut05_K2
final_DeCAF_classifier$anno$subtype[final_DeCAF_classifier$anno$subtype == "iCAF"] = "restCAF"
final_DeCAF_classifier$anno$subtype[final_DeCAF_classifier$anno$subtype == "myCAF"] = "permCAF"

## Save
save(final_DeCAF_classifier, file = "DeCAF_classifier/final_DeCAF_classifier.Rdata")



