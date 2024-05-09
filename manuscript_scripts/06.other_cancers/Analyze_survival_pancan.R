############################## Functions and libraries ##############################
# set dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load libraries
library(ConsensusClusterPlus) # R3.xx and R4.xx have different versions of ConsensusClusterPlus
library("RColorBrewer")
library(openxlsx)
library(stringr)
library(survival)
library(survminer)

# load functions
file.sources <- list.files("../R/R/",pattern="*.R")
file.sources <- paste("../R/R/", file.sources, sep="")
sapply(file.sources, source)

# load subtype info
### This is a combined subtype object with
### subtypeColList, subtypeGeneList, subtypeList and schemaList
load("../data/cmbSubtypes.RData")
print("Subtype schemas available for use:")
print(schemaList)
version_PurISS = "final"
TSPgenes <- subtypeGeneList[[33]]

# define functions
flt_sample <- function(dataSet, sampFlt) {
  dataSet$sampInfo <- dataSet$sampInfo[-sampFlt, ]
  dataSet$ex <- dataSet$ex[,-sampFlt]
  return(dataSet)
}

############################## Subtyping ################################
# list RData
dataList <- list.files("../data/TCGA_RNAseq/",pattern=".rds")

# initialize panCan gene %
#pan.mycaf <- data.frame(matrix(NA,
#                              nrow = 9,
#                              ncol = 0))
#pan.icaf <- data.frame(matrix(NA,
#                              nrow = 9,
#                              ncol = 0))

# load clinical data
clinicalAll <- read.xlsx("../data/TCGA_RNAseq/TCGA-CDR-SupplementalTableS1.xlsx")
clinicalAll$RNAseqID <- NA
clinicalAll$PurISS.prob <- NA
clinicalAll$PurISS <- NA
clinicalAll$PurISS_graded <- NA

# initial pvalue dataframe
survPvalue <- matrix(ncol = 9)
colnames(survPvalue) <- c("CancerType","N","N_events","permCAF","permCAF median","restCAF","restCAF median","HR","Pvalue")

survPvalue.DFI <- matrix(ncol = 9)
colnames(survPvalue.DFI) <- c("CancerType","N","N_events","permCAF","permCAF median","restCAF","restCAF median","HR","Pvalue")

survPvalue.DSS <- matrix(ncol = 9)
colnames(survPvalue.DSS) <- c("CancerType","N","N_events","permCAF","permCAF median","restCAF","restCAF median","HR","Pvalue")

survPvalue.PFI <- matrix(ncol = 9)
colnames(survPvalue.PFI) <- c("CancerType","N","N_events","permCAF","permCAF median","restCAF","restCAF median","HR","Pvalue")

# plot
#par(mfrow=c(2,1))
pdf(paste("../figure_DeCAF/","panCan_surv.pdf",sep=""))
for (dataName in dataList) {
  cancerType <- gsub(".rds", "", dataName)
  dataSet <- readRDS(paste("../data/TCGA_RNAseq/",cancerType,".rds",sep=""))
  sampSub <- 1:nrow(dataSet$sampInfo)
  dataSet <- Call_PurISS(dataSet, version = version_PurISS)
  
  # survival
  survDat <- dataSet$sampInfo
  tmpSplit <- data.frame(str_split_fixed(survDat$sampleID, "-", 5))
  survDat$barcode <- apply(tmpSplit[,c(1:3)], 1, paste, collapse = "-")
  survDat$loci <- tmpSplit$X4
  if(cancerType != "TCGA_LAML") {
    survDat$barcode[which(!(survDat$loci %in% "01A"))] <- paste(survDat$barcode[which(!(survDat$loci %in% "01A"))], "-",
                                                                survDat$loci[which(!(survDat$loci %in% "01A"))], sep="")
  }
  
  idxSamp <- match(survDat$barcode,clinicalAll$bcr_patient_barcode)
  survDat <- cbind(survDat,
                   clinicalAll[idxSamp,c("OS","OS.time","DFI","DFI.time","DSS","DSS.time","PFI","PFI.time")])
  survDat$OS.time <- survDat$OS.time/30
  survDat$DFI.time <- survDat$DFI.time/30
  survDat$DSS.time <- survDat$DSS.time/30
  survDat$PFI.time <- survDat$PFI.time/30
  survDat$PurISS <- factor(survDat$PurISS.final, levels = c("myCAF","iCAF"))
  
  break.time = 24
  if(cancerType == "TCGA_BRCA" | cancerType == "TCGA_BLCA" | cancerType == "TCGA_SKCM"){
    break.time = 48
  }
  # OS -------------------------------------------------------------------------
  km <- with(survDat, Surv(OS.time,OS))
  # get stats
  p = coxph(km ~ PurISS, data = survDat)
  permN <- length(which(survDat$PurISS %in% "myCAF" & !is.na(survDat$OS)))
  restN <- length(which(survDat$PurISS %in% "iCAF"  & !is.na(survDat$OS)))
  bic = round(BIC(p),3)
  hr <- paste0(round(1/summary(p)$coefficients[2],3), "(",round(1/summary(p)$conf.int[4],3), ",",round(1/summary(p)$conf.int[3],3), ")")
  pval <- round(summary(p)$logtest[3],3)

  km_fit <- survfit(km ~ PurISS, data = survDat, type = "kaplan-meier")
  permMed <- paste0(round(surv_median(km_fit)[1,2],3),"(",round(surv_median(km_fit)[1,3],3),",",round(surv_median(km_fit)[1,4],3),")")
  restMed <- paste0(round(surv_median(km_fit)[2,2],3),"(",round(surv_median(km_fit)[2,3],3),",",round(surv_median(km_fit)[2,4],3),")")
  survPvalue <- rbind(survPvalue, c(cancerType, p$n, p$nevent,permN, permMed,restN,restMed,hr,pval))
  kmplot <- ggsurvplot(km_fit, 
           conf.int = F, 
           pval = T,
           legend.title="",
           break.time.by = break.time,
           legend.labs=c("permCAF","restCAF"),
           palette = c("violetred1","turquoise4"),
           xlab = "Time (months)", 
           #ylab = "Proportion",
           pval.size = 8,
           font.title = 16,
           font.legend = 16,
           font.x = 16, 
           font.y = 16, 
           font.tickslab = 16,
           risk.table.fontsize = 7,
           censor.size=3, 
           risk.table = T,
           size = 0.3, 
           risk.table.height = 0.35,
           surv.median.line = "hv", 
           title = paste(cancerType,"\nOS HR=",hr,sep="") )
print(kmplot)

# save PurISS calls to clinical sheet
idxTmp <- which(clinicalAll$bcr_patient_barcode %in% survDat$barcode)
clinicalAll[idxTmp, c("RNAseqID","PurISS.prob.final","PurISS.final","PurISS_graded.final")] <- survDat[match(clinicalAll$bcr_patient_barcode[idxTmp], survDat$barcode),
                                                                                     c("sampleID","PurISS.prob.final","PurISS.final","PurISS_graded.final")]
clinicalAll$RNAseqID[idxTmp] <- as.character(survDat$sampleID[match(clinicalAll$bcr_patient_barcode[idxTmp], survDat$barcode)])

# DFI -------------------------------------------------------------------------
km <- with(survDat, Surv(DFI.time,DFI))
# get stats
if(cancerType == "TCGA_LAML" | cancerType == "TCGA_GBM" | cancerType == "TCGA_SKCM"| cancerType == "TCGA_THYM"| cancerType == "TCGA_UVM") {
  survPvalue.DFI <- rbind(survPvalue.DFI, c(cancerType, nrow(survDat), NA,permN, NA,restN,NA,NA, NA))
} else {
p = coxph(km ~ PurISS, data = survDat)
permN <- length(which(survDat$PurISS %in% "myCAF" & !is.na(survDat$DFI)))
restN <- length(which(survDat$PurISS %in% "iCAF"  & !is.na(survDat$DFI)))
bic = round(BIC(p),3)
hr <- paste0(round(1/summary(p)$coefficients[2],3), "(",round(1/summary(p)$conf.int[4],3), ",",round(1/summary(p)$conf.int[3],3), ")")
pval <- round(summary(p)$logtest[3],3)
km_fit <- survfit(km ~ PurISS, data = survDat, type = "kaplan-meier")
permMed <- paste0(round(surv_median(km_fit)[1,2],3),"(",round(surv_median(km_fit)[1,3],3),",",round(surv_median(km_fit)[1,4],3),")")
restMed <- paste0(round(surv_median(km_fit)[2,2],3),"(",round(surv_median(km_fit)[2,3],3),",",round(surv_median(km_fit)[2,4],3),")")
survPvalue.DFI <- rbind(survPvalue.DFI, c(cancerType, p$n, p$nevent,permN, permMed,restN,restMed,hr,pval))
kmplot <- ggsurvplot(km_fit, 
                     conf.int = F, 
                     pval = T,
                     #pval.coord = c(0.5*max(survDat$OS.time[which(!is.na(survDat$OS.time))]), 0.9),
                     legend.title="",
                     break.time.by = break.time,
                     legend.labs=c("permCAF","restCAF"),
                     palette = c("violetred1","turquoise4"),
                     xlab = "Time (months)", 
                     #ylab = "Proportion",
                     pval.size = 8,
                     font.title = 16,
                     font.legend = 16,
                     font.x = 16, 
                     font.y = 16, 
                     font.tickslab = 16,
                     risk.table.fontsize = 7,
                     censor.size=3, 
                     risk.table = T,
                     size = 0.3, 
                     risk.table.height = 0.35,
                     surv.median.line = "hv", 
                     title = paste(cancerType,"\nDFI HR=",hr,sep="") )
print(kmplot)
}

# DSS -------------------------------------------------------------------------
km <- with(survDat, Surv(DSS.time,DSS))
# get stats
if(cancerType == "TCGA_LAML") {
  survPvalue.DSS <- rbind(survPvalue.DSS, c(cancerType, nrow(survDat), NA,permN, NA,restN,NA,NA, NA))
} else {
  p = coxph(km ~ PurISS, data = survDat)
  permN <- length(which(survDat$PurISS %in% "myCAF" & !is.na(survDat$DSS)))
  restN <- length(which(survDat$PurISS %in% "iCAF"  & !is.na(survDat$DSS)))
  bic = round(BIC(p),3)
  hr <- paste0(round(1/summary(p)$coefficients[2],3), "(",round(1/summary(p)$conf.int[4],3), ",",round(1/summary(p)$conf.int[3],3), ")")
  pval <- round(summary(p)$logtest[3],3)
  km_fit <- survfit(km ~ PurISS, data = survDat, type = "kaplan-meier")
  permMed <- paste0(round(surv_median(km_fit)[1,2],3),"(",round(surv_median(km_fit)[1,3],3),",",round(surv_median(km_fit)[1,4],3),")")
  restMed <- paste0(round(surv_median(km_fit)[2,2],3),"(",round(surv_median(km_fit)[2,3],3),",",round(surv_median(km_fit)[2,4],3),")")
  survPvalue.DSS <- rbind(survPvalue.DSS, c(cancerType, p$n, p$nevent,permN, permMed,restN,restMed,hr,pval))
  kmplot <- ggsurvplot(km_fit, 
                       conf.int = F, 
                       pval = T,
                       #pval.coord = c(0.5*max(survDat$OS.time[which(!is.na(survDat$OS.time))]), 0.9),
                       legend.title="",
                       break.time.by = break.time,
                       legend.labs=c("permCAF","restCAF"),
                       palette = c("violetred1","turquoise4"),
                       xlab = "Time (months)", 
                       #ylab = "Proportion",
                       pval.size = 8,
                       font.title = 16,
                       font.legend = 16,
                       font.x = 16, 
                       font.y = 16, 
                       font.tickslab = 16,
                       risk.table.fontsize = 7,
                       censor.size=3, 
                       risk.table = T,
                       size = 0.3, 
                       risk.table.height = 0.35,
                       surv.median.line = "hv", 
                       title = paste(cancerType,"\nDSS HR=",hr,sep="") )
  print(kmplot)
}

# PFI -------------------------------------------------------------------------
km <- with(survDat, Surv(PFI.time,PFI))
# get stats
if(cancerType == "TCGA_LAML") {
  survPvalue.PFI <- rbind(survPvalue.PFI, c(cancerType, nrow(survDat), NA,permN, NA,restN,NA,NA, NA))
} else {
  p = coxph(km ~ PurISS, data = survDat)
  permN <- length(which(survDat$PurISS %in% "myCAF" & !is.na(survDat$PFI)))
  restN <- length(which(survDat$PurISS %in% "iCAF"  & !is.na(survDat$PFI)))
  bic = round(BIC(p),3)
  hr <- paste0(round(1/summary(p)$coefficients[2],3), "(",round(1/summary(p)$conf.int[4],3), ",",round(1/summary(p)$conf.int[3],3), ")")
  pval <- round(summary(p)$logtest[3],3)
  km_fit <- survfit(km ~ PurISS, data = survDat, type = "kaplan-meier")
  permMed <- paste0(round(surv_median(km_fit)[1,2],3),"(",round(surv_median(km_fit)[1,3],3),",",round(surv_median(km_fit)[1,4],3),")")
  restMed <- paste0(round(surv_median(km_fit)[2,2],3),"(",round(surv_median(km_fit)[2,3],3),",",round(surv_median(km_fit)[2,4],3),")")
  survPvalue.PFI <- rbind(survPvalue.PFI, c(cancerType, p$n, p$nevent,permN, permMed,restN,restMed,hr,pval))
  kmplot <- ggsurvplot(km_fit, 
                       conf.int = F, 
                       pval = T,
                       #pval.coord = c(0.5*max(survDat$OS.time[which(!is.na(survDat$OS.time))]), 0.9),
                       legend.title="",
                       break.time.by = break.time,
                       legend.labs=c("permCAF","restCAF"),
                       palette = c("violetred1","turquoise4"),
                       xlab = "Time (months)", 
                       #ylab = "Proportion",
                       pval.size = 8,
                       font.title = 16,
                       font.legend = 16,
                       font.x = 16, 
                       font.y = 16, 
                       font.tickslab = 16,
                       risk.table.fontsize = 7,
                       censor.size=3, 
                       risk.table = T,
                       size = 0.3, 
                       risk.table.height = 0.35,
                       surv.median.line = "hv", 
                       title = paste(cancerType,"\nPFI HR=",hr,sep="") )
  print(kmplot)
}
}
dev.off()

# 
survPvalue <- as.data.frame(survPvalue[-1,])
survPvalue.DFI <- as.data.frame(survPvalue.DFI[-1,])
survPvalue.DSS <- as.data.frame(survPvalue.DSS[-1,])
survPvalue.PFI <- as.data.frame(survPvalue.PFI[-1,])

survPvalue$Pvalue.bonferroni <- p.adjust(survPvalue$Pvalue, method = "bonferroni")
survPvalue$Pvalue.fdr <- p.adjust(survPvalue$Pvalue, method = "fdr")
survPvalue$Pvalue.bh <- p.adjust(survPvalue$Pvalue, method = "BH")
write.xlsx(survPvalue, "../figure_DeCAF/pancan_os_pvalue.xlsx")

survPvalue.DFI$Pvalue.bonferroni <- p.adjust(survPvalue.DFI$Pvalue, method = "bonferroni")
survPvalue.DFI$Pvalue.fdr <- p.adjust(survPvalue.DFI$Pvalue, method = "fdr")
survPvalue.DFI$Pvalue.bh <- p.adjust(survPvalue.DFI$Pvalue, method = "BH")
write.xlsx(survPvalue.DFI, "../figure_DeCAF/pancan_dfi_pvalue.xlsx")

survPvalue.DSS$Pvalue.bonferroni <- p.adjust(survPvalue.DSS$Pvalue, method = "bonferroni")
survPvalue.DSS$Pvalue.fdr <- p.adjust(survPvalue.DSS$Pvalue, method = "fdr")
survPvalue.DSS$Pvalue.bh <- p.adjust(survPvalue.DSS$Pvalue, method = "BH")
write.xlsx(survPvalue.DSS, "../figure_DeCAF/pancan_dss_pvalue.xlsx")

survPvalue.PFI$Pvalue.bonferroni <- p.adjust(survPvalue.PFI$Pvalue, method = "bonferroni")
survPvalue.PFI$Pvalue.fdr <- p.adjust(survPvalue.PFI$Pvalue, method = "fdr")
survPvalue.PFI$Pvalue.bh <- p.adjust(survPvalue.PFI$Pvalue, method = "BH")
write.xlsx(survPvalue.PFI, "../figure_DeCAF/pancan_pfi_pvalue.xlsx")

names(clinicalAll)[match(c("PurISS.prob","PurISS","PurISS_graded"),names(clinicalAll))] <- c("DeCAF_prob","DeCAF","DeCAF_graded")
write.xlsx(clinicalAll, "../figure_DeCAF/pancan_clinical_data_decaf.xlsx")
