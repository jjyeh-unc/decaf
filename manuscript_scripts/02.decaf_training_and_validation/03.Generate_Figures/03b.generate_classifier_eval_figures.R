#####################################################################################

##########     Classifier Evaluation Figure Generating Script             ############

####################################################################################




######################################################################################
############## Step 1: Create Functions for Generating Figures ########################
######################################################################################




# load libraries
library(pheatmap)
library(ggplot2)
library(scales)
library(caret)
library(pROC)
library(kableExtra)
library(nsROC)
library(ROCR)

## Colors for DeCAF
DeCAFList = c("permCAF","restCAF")
DeCAFCol = c("violetred1","turquoise4")

############################################################

## Function for Heatmap
Plot_heatmap_function <- function(subtype_data, gene_data, DeCAF_data, classifier, training, CAF_type = c("K2", "K4"), schemaList, include_PurIST = F){
  
  ## List of Samples To Include
  if(CAF_type == "K4"){
    samples_subset = !(subtype_data$dataSet %in% training) & !(is.na(subtype_data$SCISSORS_CAF_ALL))
  }
  
  if(CAF_type == "K2"){
    samples_subset = !(subtype_data$dataSet %in% training) & !(is.na(subtype_data$SCISSORS_CAF))
  }
  
  
  ## Transform TSP's to ranks
  TSPgeneMat = apply(gene_data[,samples_subset], 2, rank)
  
  ## Transfrom Calls to Ordered Factor
  DeCAF_data$Stroma_Subtype_graded = factor(DeCAF_data$Stroma_Subtype_graded, levels = c("Strong restCAF", "Likely restCAF", "Lean restCAF",
                                                                                           "Lean permCAF", "Likely permCAF", "Strong permCAF"))
  ## Make ColorScapr
  colscale = colorRampPalette(c("turquoise4","white", "violetred1"))(6)
  colscale2 = colorRampPalette(c("turquoise4","white", "violetred1"))(5)
  colscale3 = colorRampPalette(c("turquoise4","white", "violetred1"))(2)
  colscale4 = c(colscale2, "black")
  
  # Create cleaned up study labels
  Study = subtype_data$dataSet[samples_subset]
  Study = gsub("_Pi","",Study)
  
  
  
  ## Change CAF based on if K2 or K4 required
  if(CAF_type == "K4"){
    cc_subtype = subtype_data$SCISSORS_CAF_ALL[samples_subset]
    cc_subtype_names = c("restCAF", "Mixed.restCAF","Mixed","Mixed.permCAF", "permCAF", "Absent")
    cc_colscale = colscale4
  }
  if(CAF_type == "K2"){
    cc_subtype = subtype_data$SCISSORS_CAF[samples_subset]
    cc_subtype_names = c("restCAF", "permCAF")
    cc_colscale = colscale3
  }
  
  
  # create annotation vector for the samples
  annotation_col = data.frame(Study = Study, DeCAF = DeCAF_data$Stroma_Subtype_graded[samples_subset],
                              'CC Subtype' = cc_subtype,
                              # "Elyada CAF" = subtype_data$Elyada_CAF[samples_subset],
                              #"Moffitt Stroma" = subtype_data$MS_K2[samples_subset], 
                              check.names = F)
  if(include_PurIST == T){
    annotation_col$PurIST = subtype_data$PurIST[samples_subset]
  }
  rownames(annotation_col) = colnames(gene_data)[samples_subset]
  
  # create annotation color vector

  index_PurIST = which(schemaList == "PurIST")
  annotation_col_colors = list(
    DeCAF = colscale,
    'CC Subtype' = cc_colscale,
    'Training Labels' = DeCAFCol,
    'Up Subtype' =DeCAFCol
  )
  names(annotation_col_colors$'CC Subtype')  = cc_subtype_names
  names(annotation_col_colors$'Training Labels') = DeCAFList
  names(annotation_col_colors$DeCAF) = levels(DeCAF_data$Stroma_Subtype_graded)
  names(annotation_col_colors$'Up Subtype') = DeCAFList
  
  if(include_PurIST == T){
    annotation_col_colors$PurIST = subtypeColList[[index_PurIST]]
    names(annotation_col_colors$PurIST) = subtypeList[[index_PurIST]]
  }
  
  
  
  ## Sort columns and rows of heatmap
  sort = order(DeCAF_data$Pred_prob_permCAF[samples_subset])
  tsp = match(as.vector(t(classifier$TSPs[classifier$fit$beta[-1]!=0,])),rownames(TSPgeneMat))
  
  ## Classifier gene annotation
  annotation_row = data.frame('Up Subtype' = factor((rep(DeCAFList  ,length(tsp)/2))),check.names = F)
  rownames(annotation_row) = rownames(TSPgeneMat)[tsp]
  
  
  ## Generate heatmap
  heatmap = pheatmap(TSPgeneMat[tsp,sort], annotation_col = annotation_col[sort,],annotation_row = annotation_row, 
                     annotation_colors = annotation_col_colors,show_colnames = F, cluster_cols = F,cluster_rows = F,
                     annotation_names_row = F,width = 12, height = 12,fontsize = 11)
  
  
  
  ## Output
  return(heatmap)
  
  
}

#########################################################################

## Function for Barplot of Strong vs Likely vs Lean
# 4C Generating Function


Plot_graded_barplot <- function(subtype_data, DeCAF_data, classifier, training){
  
  

## Subset of Samples to Remove
samples_subset = !(subtype_data$dataSet %in% training) & subtype_data$CAF_whitelist & subtype_data$dataSet != "YehSeq_Pi"


## Relevel DeCAF Graded
DeCAF_data$Stroma_Subtype_graded = factor(DeCAF_data$Stroma_Subtype_graded, levels = c("Strong restCAF", "Likely restCAF", "Lean restCAF",
                                                                                         "Lean permCAF", "Likely permCAF", "Strong permCAF"))
## Make ColorScale
colscale =  colorRampPalette(c("turquoise4","white", "violetred1"))(6)


## Make barplot
Caf_pred = DeCAF_data$Stroma_Subtype_graded[samples_subset]
Caf_pred = data.frame(Caf_pred)
Caf_pred2 = data.frame(DeCAF_prob = DeCAF_data$Pred_prob_permCAF[samples_subset], Scissors_CAF = subtype_data$SCISSORS_CAF_ALL[samples_subset])
library(ggplot2)
colscale2 = colorRampPalette(c("turquoise4","white", "violetred1"))(5)
Caf_pred2 = Caf_pred2[Caf_pred2$Scissors_CAF != "Absent",]
Caf_pred2$Scissors_CAF = factor(Caf_pred2$Scissors_CAF , levels = c("restCAF", "Mixed.restCAF","Mixed","Mixed.permCAF", "permCAF"))

barplot = ggplot(Caf_pred, aes(x = Caf_pred, fill = Caf_pred)) +geom_bar(aes(y=after_stat(count)/sum(after_stat(count))))+ 
  theme_classic() + scale_fill_manual(values = colscale) +
  scale_y_continuous(labels =scales::percent, limits = c(0,0.5))+ labs(y = "Percentage of Samples", x = "")+ 
  theme(legend.position = "none", axis.text.y = element_text(color= "black"),axis.text.x = element_text(angle = 45, vjust = 0.5, color = "black"), text = element_text(size=22, color = "black")) 
barplot

return(barplot)
}

#############################################################################

## Function for barplot of lean likely strong CAF classification with proportion of CC Calls 
Plot_graded_Proportion_barplot <- function(subtype_data, DeCAF_data, classifier, training, CAF_type = c("K2", "K4")){
  
  
  ## Subset of Samples to Remove
  samples_subset = !(subtype_data$dataSet %in% training) & subtype_data$CAF_whitelist
  
  ## Relevel DeCAF Graded
  DeCAF_data$Stroma_Subtype_graded = factor(DeCAF_data$Stroma_Subtype_graded, levels = c("Strong restCAF", "Likely restCAF", "Lean restCAF",
                                                                                           "Lean permCAF", "Likely permCAF", "Strong permCAF"))
  ## Make ColorScale
  colscale = colorRampPalette(c("turquoise4","white", "violetred1"))(6)
  
  ## Only Keep Subset Samples
  DeCAF_data = DeCAF_data[samples_subset,]
  subtype_data = subtype_data[samples_subset,]
  
  ## Obtain CAF Col, Colors, and Names
  if(CAF_type == "K2"){
    CC_calls = subtype_data$SCISSORS_CAF
    cc_subtype_names = c("restCAF", "permCAF")
    cc_colscale = colorRampPalette(c("turquoise4","white", "violetred1"))(2)
  }
  if(CAF_type == "K4"){
    CC_calls = subtype_data$SCISSORS_CAF_ALL
    CC_calls = factor(CC_calls, level = c("restCAF", "Mixed.restCAF","Mixed","Mixed.permCAF", "permCAF", "Absent"))
    cc_subtype_names = c("restCAF", "Mixed.restCAF","Mixed","Mixed.permCAF", "permCAF", "Absent")
    cc_colscale = colorRampPalette(c("turquoise4","white", "violetred1"))(5)
    cc_colscale[6] = "black"
  }
  
  
  
  

  ## Tabulate
  tab = table(CC_calls,DeCAF_data$Stroma_Subtype_graded)
  values = c(tab)
  CC_calls_df = rep(cc_subtype_names, length(levels(DeCAF_data$Stroma_Subtype_graded)))
  DeCAF_calls_df = rep(levels(DeCAF_data$Stroma_Subtype_graded), each = length(cc_subtype_names))
  DeCAF_calls_df = factor(DeCAF_calls_df, levels = levels(DeCAF_data$Stroma_Subtype_graded))
  df = data.frame(values = values, CC_calls = CC_calls_df, DeCAF_calls = DeCAF_calls_df)
  
  if(CAF_type == "K4"){
    df$CC_calls = factor(df$CC_calls, level = c("restCAF", "Mixed.restCAF","Mixed","Mixed.permCAF", "permCAF", "Absent"))
      }
  
  
  ## Make Barplot
  barplot = ggplot(df, aes(fill = CC_calls, y = values, x = DeCAF_calls))+
    geom_bar(position = "fill", stat = "identity", color = "black") + theme_classic()  +
    scale_fill_manual(values = cc_colscale, labels = cc_subtype_names, name = "") +
    labs(y = "Subtype Percent", x = "")+ 
    theme(legend.position = "bottom", axis.text.y = element_text(color = "black"), axis.text.x = element_text(angle = 45, vjust = 0.5, color = "black"), text = element_text(size=22))+
    scale_y_continuous(labels =scales::percent)
  
  return(barplot)
  
}

###################################################################################################

## Function to Generate Pooled ROC
Plot_Pooled_ROC <- function(subtype_data, DeCAF_data,  training){
  
  ## Testing Set
  testing0 =  unique(subtype_data$dataSet)[!unique(subtype_data$dataSet) %in% training]
  
  ## Sample Subset
  samples_subset = !(subtype_data$dataSet %in% training) & !is.na(subtype_data$SCISSORS_CAF)
  DeCAF_data = DeCAF_data[samples_subset,]
  subtype_data = subtype_data[samples_subset,]
  
  ## Calculate TP, FP, TN, FN
  data = NULL
  for(i in 1:(length(testing0))){
    pred = DeCAF_data$Stroma_Subtype[DeCAF_data$dataSet == testing0[i]]
    pred_p = DeCAF_data$Pred_prob_permCAF[DeCAF_data$dataSet == testing0[i]]
    label = subtype_data$SCISSORS_CAF[subtype_data$dataSet == testing0[i]]
    pred = prediction(prediction = pred_p, labels = as.numeric(label == "permCAF"))
    TP = as.numeric(attr(pred,"tp")[[1]])
    FP = as.numeric(attr(pred,"fp")[[1]])
    TN = as.numeric(attr(pred,"tn")[[1]])
    FN = as.numeric(attr(pred,"fn")[[1]])
    Author = rep(i, length(TP))
    data = rbind(data, cbind(Author, TP, FP, TN, FN))
  }
  data  = data.frame(data)
  
  
  
  
  
    ## Make ROC curve
    return(
      metaROC(data,cex.main = 1,cex.lab =1,cex.axis = 1, model="random-effects", plot.Author=TRUE,plot.inter.var = T,cex.Author = .5,lwd.Author = 2,alpha.trans = .25,Ni = 100)
      
    )  
  
  
  
}
  


#####################################################################################################


## Function to Generate Table of Classifier Performance
Table_Classif_Perform <- function(subtype_data, DeCAF_data,  training){
  
  ## Testing Set
  testing0 =  unique(subtype_data$dataSet)[!unique(subtype_data$dataSet) %in% training]
  
  ## Sample Subset
  samples_subset = !(subtype_data$dataSet %in% training) & !is.na(subtype_data$SCISSORS_CAF)
  DeCAF_data = DeCAF_data[samples_subset,]
  subtype_data = subtype_data[samples_subset,]
  
  ## Empty Vector for metrics
  N = c()
  permCAF = c()
  Accuracy = c()
  Sensitivity = c()
  Specificity = c()
  AUC = c()
  
  
  samples_subset_train = (SubtypeCmb$dataSet %in% training) & !is.na(SubtypeCmb$SCISSORS_CAF)
  DeCAF_data_t = DeCAFCmb[samples_subset_train,]
  subtype_data_t = SubtypeCmb[samples_subset_train,]
  
  
  ## Obtain Metric for each study 
  for(i in 1:(length(testing0))){
    pred = DeCAF_data$Stroma_Subtype[DeCAF_data$dataSet == testing0[i]]
    pred_p = DeCAF_data$Pred_prob_permCAF[DeCAF_data$dataSet == testing0[i]]
    label = subtype_data$SCISSORS_CAF[subtype_data$dataSet == testing0[i]]
    N[i] = length(label)
    permCAF[i] = sum(label == "permCAF")
    Accuracy[i] = round(mean(pred == label),3)
    Sensitivity[i] = round(sensitivity(as.factor(pred), as.factor(label), positive = "permCAF"),3)
    Specificity[i] = round(specificity(as.factor(pred), as.factor(label), negative = "restCAF"),3)
    roc = roc(response = label, predictor = pred_p, levels = c("restCAF", "permCAF"), quiet = T)
    AUC[i] = round(roc$auc[1],3)
    
    
  }
  
  ## Combien Results into one df
  pmat = data.frame(N, permCAF, Accuracy, Sensitivity, Specificity, AUC)
  
  ## Clean up Study Name
  testing0 = gsub("_PrimaryPDAC","",testing0)
  rownames(pmat) = testing0[1:(length(testing0))]
  
  
  
  return(pmat)
  

  
}



  
######################################################################################
##############          Step 2: Load In Data                 ########################
######################################################################################

## Set working directory
## setwd(PATH TO DeCAF_Classifier_Clean folder)

## Load Subtype, Gene Data, and DeCAF Calls
SubtypeCmb <- read.csv("data/clean_data/combined_data/combined_subtype_DeCAF.csv")
GeneDataCmb <- read.csv("data/clean_data/combined_data/combined_classifier_gene_exp_DeCAF.csv",
                        row.names = 1)
DeCAFCmb <- read.csv("data/clean_data/combined_data/combined_classifier_calls_DeCAF.csv")

## Load Classifier to obtain training groups
load("DeCAF_classifier/final_DeCAF_classifier.Rdata")
training = unique(final_DeCAF_classifier$anno$study)
classifier = final_DeCAF_classifier$classifier2


## Load Schema List - Aka Color Codes for other classifiers
load("data/color_codes/cmbsubtypes.RData")



######################################################################################
##############          Step 3: Generate Figures              ########################
######################################################################################

#CAF K4 Figure Generator
K4_heatmap = Plot_heatmap_function(subtype_data = SubtypeCmb, 
                                   gene_data = GeneDataCmb, 
                                   DeCAF_data = DeCAFCmb, 
                                   classifier = classifier, 
                                   training = training, 
                                   CAF_type = "K4",
                                   schemaList = schemaList, include_PurIST = T)

filename =  paste("data/results/figures/DeCAF_model_figures/ALL_w_PurIST.pdf", sep ="")
pdf(filename, width = 10, height = 7)
K4_heatmap
dev.off()


#CAF K2 Figure Generator
K2_heatmap = Plot_heatmap_function(subtype_data = SubtypeCmb, 
                                   gene_data = GeneDataCmb, 
                                   DeCAF_data = DeCAFCmb, 
                                   classifier = classifier, 
                                   training = training, 
                                   CAF_type = "K2",
                                   schemaList = schemaList, include_PurIST = T)





## Figure Graded Barplot
grad_barplot <- Plot_graded_barplot(subtype_data = SubtypeCmb, DeCAF_data = DeCAFCmb, classifier  = classifier, training = training)
filename =  paste("data/results/figures/DeCAF_model_figures/Graded_Barplot.pdf", sep ="")
pdf(filename)
grad_barplot
dev.off()



## Figure Proportion Barplot
Proportion_Barplot <- Plot_graded_Proportion_barplot(subtype_data = SubtypeCmb, DeCAF_data = DeCAFCmb, classifier  = classifier, training = training, 
                                           CAF_type = "K4")

filename =  paste("data/results/figures/DeCAF_model_figures/Proportion_Barplot.pdf", sep ="")
pdf(filename)
Proportion_Barplot
dev.off()


## Figure 4E Generator
ROC_curve = Plot_Pooled_ROC(subtype_data = SubtypeCmb, DeCAF_data = DeCAFCmb,  training = training)
filename =  paste("data/results/figures/DeCAF_model_figures/Meta_ROC.pdf", sep ="")

pdf(filename, height = 7, width = 7, pointsize = 25)
Plot_Pooled_ROC(subtype_data = SubtypeCmb, DeCAF_data = DeCAFCmb,  training = training)
dev.off()
dev.off()
## Click and Save Pop Up Window!
ROC_curve = Plot_Pooled_ROC(subtype_data = SubtypeCmb, DeCAF_data = DeCAFCmb,  training = training)


## Figure 4G Generator
table = Table_Classif_Perform(subtype_data = SubtypeCmb, DeCAF_data = DeCAFCmb,  training = training)
write.csv(table, "data/results/tables/Metrics_ByStudy.csv")











