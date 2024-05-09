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
library(reshape2)
library(dplyr)
library(patchwork)
# library(forcats) fct_rev()

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

# load dataset
Load_cafSubtype <- function(rDataName) {
  survDat <- readRDS(paste("../data_asis_refseq_star_rsem/",rDataName,".caf_subtype.rds",sep=""))
  return(survDat)
}

Load_survDat <- function(rDataName) {
  survDat <- readRDS(paste("../data_asis_refseq_star_rsem/",rDataName,".survival_data.rds",sep=""))
  return(survDat)
}

# define fuctions
plot_surv <- function(survDat, km, cancerType) {
  # get stats
  p = coxph(km ~ DeCAF, data = survDat)
  bic = round(BIC(p),3)
  hr <- round(1/summary(p)$coefficients[2],3)
  ci_upper <- round(1/summary(p)$conf.int[3],3)
  ci_lower <- round(1/summary(p)$conf.int[4],3)
  pval <- round(summary(p)$sctest[3],3)
  survPvalue <- c(cancerType, 
                  p$n, p$nevent,
                  hr, ci_lower, ci_upper,
                  pval)
  
  # plot km
  km_fit <- survfit(km ~ DeCAF, data = survDat, type = "kaplan-meier")
  kmplot <- ggsurvplot(font.legend = 16, 
                       font.title = 20,        
                       font.x = 16,  
                       font.y = 16, 
                       font.tickslab = 14,
                       risk.table.fontsize = 6,
                       risk.table.height = 0.35,
                       surv.median.line = "hv", 
                       size = 0.3, 
                       censor.size=7, 
                       pval.size = 8,
                       km_fit, 
                       conf.int = F, 
                       pval = T,
                       pval.coord = c(0.5*max(survDat$time[which(!is.na(survDat$time))]), 0.75),
                       risk.table = T,
                       legend.title="",
                       break.time.by = 12,
                       legend.labs=c("permCAF","restCAT"),
                       palette = c("violetred1","turquoise4"),
                       xlab = "Time (months)", 
                       title = paste(gsub("TCGA_","",survPvalue[1]),
                                     " HR=",survPvalue[4],
                                     " (",survPvalue[5],",",survPvalue[6],")",sep="") )
  print(kmplot)
}

############################## Subtyping ################################
#clinicalAll <- read.xlsx("../data/TCGA_RNAseq/TCGA-CDR-SupplementalTableS1.xlsx")

survPvalue <- read.xlsx("../results/caper_pancan_os_pvalue.xlsx")

# update TCGA PAAD to N=146
rDataName <- "TCGA_PAAD"
cafSubtype <- Load_cafSubtype(rDataName)

# call PurISS
dataSet <- readRDS(paste("../data_asis_refseq_star_rsem/",rDataName,".rds",sep=""))
## PurISS
version_PurISS = "final"
dataSet <- Call_PurISS(dataSet, version = version_PurISS)
cafSubtype$Subtype[paste("PurISS", sep = ".")] <- dataSet$sampInfo[[paste("PurISS", version_PurISS, sep = ".")]]
cafSubtype$Subtype[paste("PurISS_graded", sep = ".")] <- dataSet$sampInfo[[paste("PurISS_graded", version_PurISS, sep = ".")]]
cafSubtype$Subtype[paste("PurISS.prob", sep = ".")] <- dataSet$sampInfo[[paste("PurISS.prob", version_PurISS, sep = ".")]]

# parse cafSubtype and survDat
cafSubtype <- cafSubtype$Subtype
survDat <- Load_survDat(rDataName)
survDat <- survDat[,c("sampID","time","event","whitelist")]
survDat$time <- as.numeric(survDat$time)
survDat$event <- as.numeric(survDat$event)
survDat$Dataset <- rDataName
survDat <- cbind(survDat, cafSubtype)
survDat <- survDat[which(survDat$whitelist),]
survDat$PurISS <- factor(survDat$PurISS, levels = c("myCAF","iCAF"))

# relabel
survDat$DeCAF <- as.character(survDat$PurISS)
survDat$DeCAF[which(survDat$DeCAF %in% "myCAF")] <- "permCAF"
survDat$DeCAF[which(survDat$DeCAF %in% "iCAF")] <- "restCAF"
survDat$DeCAF <- factor(survDat$DeCAF, levels = c("permCAF","restCAF"))

# km
km <- with(survDat, Surv(time,event))

# get stats
p = coxph(km ~ PurISS, data = survDat)
permN <- length(which(survDat$DeCAF %in% "permCAF"))
restN <- length(which(survDat$DeCAF %in% "restCAF"))
bic = round(BIC(p),3)
hr <- round(1/summary(p)$coefficients[2],3)
ci_lower <- round(1/summary(p)$conf.int[4],3)
ci_upper <- round(1/summary(p)$conf.int[3],3)
pval <- round(summary(p)$sctest[3],3)
survPvalue[which(survPvalue$CancerType %in% c("TCGA_PAAD")),2:7] <- c(p$n, p$nevent,permN, restN, hr, ci_lower, ci_upper, pval)

# correct
survPvalue$Pvalue.bonferroni <- p.adjust(survPvalue$Pvalue, method = "bonferroni")
survPvalue$Pvalue.fdr <- p.adjust(survPvalue$Pvalue, method = "fdr")
survPvalue$Pvalue.bh <- p.adjust(survPvalue$Pvalue, method = "BH")
#write.xlsx(survPvalue, "../results/caper_pancan_os_pvalue.paad_revised.xlsx")

# survival
pdf("../figure_DeCAF/surv_kirc_meso_thca.pdf")
plot_surv(survDat, km, 
          cancerType = rDataName)

# plot others
clinicalAll <- read.xlsx("../data/TCGA_RNAseq/TCGA-CDR-SupplementalTableS1.xlsx")
clinicalAll$RNAseqID <- NA

for (cancerType in c("TCGA_KIRC","TCGA_MESO","TCGA_THCA")) {
  dataSet <- readRDS(paste("../data/TCGA_RNAseq/",cancerType,".rds",sep=""))
  sampSub <- 1:nrow(dataSet$sampInfo)
  dataSet <- Call_PurISS(dataSet, version = version_PurISS)
  
  # survival
  survDat <- dataSet$sampInfo
  tmpSplit <- data.frame(str_split_fixed(survDat$sampleID, "-", 5))
  survDat$barcode <- apply(tmpSplit[,c(1:3)], 1, paste, collapse = "-")
  survDat$barcode0 <- survDat$barcode
  survDat$loci <- tmpSplit$X4
  
  # 01 is solide primary tumor
  #survDat$barcode[which(!(survDat$loci %in% "01A"))] <- paste(survDat$barcode[which(!(survDat$loci %in% "01A"))], "-",
  #                                                            survDat$loci[which(!(survDat$loci %in% "01A"))], sep="")
  idxAll <- 1:nrow(survDat)
  idxGood <- grep("01",survDat$loci)
  idxBad <- idxAll[which(!(idxAll %in% idxGood))]
  survDat$barcode[idxBad] <- paste(survDat$barcode[idxBad], "-", survDat$loci[idxBad], sep="")
  
  idxSamp <- match(survDat$barcode,clinicalAll$bcr_patient_barcode)
  survDat <- cbind(survDat,
                   clinicalAll[idxSamp,c("OS","OS.time","PFI","PFI.time")])
  survDat$OS.time <- survDat$OS.time/30
  survDat$PFI.time <- survDat$PFI.time/30

  #survDat$Dataset <- cancerType
  survDat$DeCAF <- survDat$PurISS.final
  survDat$DeCAF[which(survDat$DeCAF %in% "myCAF")] <- "permCAF"
  survDat$DeCAF[which(survDat$DeCAF %in% "iCAF")] <- "restCAF"
  survDat$DeCAF <- factor(survDat$DeCAF, levels = c("permCAF","restCAF"))
  survDat$DeCAF_graded <- survDat$PurISS_graded.final
  survDat$DeCAF_graded <- gsub("myCAF","permCAF",survDat$DeCAF_graded)
  survDat$DeCAF_graded <- gsub("iCAF","restCAF",survDat$DeCAF_graded)
  survDat$DeCAF.prob <- survDat$PurISS.prob.final
  
  # print KIRC 
  #survDat$Histological.grade <- clinicalAll$histological_grade[match(survDat$barcode0,clinicalAll$bcr_patient_barcode)]
  #survDat$Survival_analysis_whitelist <- FALSE
  #survDat$Survival_analysis_whitelist[which(!is.na(survDat$OS))] <- TRUE
  #survDat$OS <- clinicalAll$OS[match(survDat$barcode0,clinicalAll$bcr_patient_barcode)]
  #survDat$OS.time <- clinicalAll$OS.time[match(survDat$barcode0,clinicalAll$bcr_patient_barcode)]
  #write.xlsx(survDat[,c("Dataset","sampleID","DeCAF" ,"DeCAF_graded","DeCAF.prob","Histological.grade","Survival_analysis_whitelist","OS","OS.time")],
  #           file = paste0("../figure_DeCAF/suppl_",cancerType,".xlsx"))
  
  # print MESO
  #survDat$Histological.type <- clinicalAll$histological_type[match(survDat$barcode0,clinicalAll$bcr_patient_barcode)]
  #survDat$Survival_analysis_whitelist <- FALSE
  #survDat$Survival_analysis_whitelist[which(!is.na(survDat$OS.time))] <- TRUE
  #survDat$OS <- clinicalAll$OS[match(survDat$barcode0,clinicalAll$bcr_patient_barcode)]
  #survDat$OS.time <- clinicalAll$OS.time[match(survDat$barcode0,clinicalAll$bcr_patient_barcode)]
  # match pathCall
  #pathCall <- read.xlsx("../data/TCGA_pancan_clinical_caf11523_AI_231114.xlsx")
  #idxSamp <- match(survDat$barcode0, pathCall$bcr_patient_barcode)
  #survDat$`Myxoid(yes.1,.no.0)` <- pathCall$`Myxoid(yes.1,.no.0)`[idxSamp]
  #survDat$`Fibrous(yes.1,.no.0)` <- pathCall$`Collagenized(yes.1,.no.0)`[idxSamp]
  #survDat$`Myx1/Mixed2/Fib3` <- pathCall$`Myx1/mixed2/Coll3`[idxSamp]
  #write.xlsx(survDat[,c("Dataset","sampleID","DeCAF" ,"DeCAF_graded","DeCAF.prob","Histological.type","Myxoid(yes.1,.no.0)","Fibrous(yes.1,.no.0)","Myx1/Mixed2/Fib3","Survival_analysis_whitelist","OS","OS.time")],
  #           file = paste0("../figure_DeCAF/suppl_",cancerType,".xlsx"))
  #write.xlsx(survDat[,c(1,12:15,17:ncol(survDat),16,8,9)],
  #           file = paste0("../figure_DeCAF/suppl_",cancerType,".xlsx"))
  
  # km
  survDat$time <- survDat$OS.time
  survDat$event <- survDat$OS
  km <- with(survDat, Surv(time,event))
  plot_surv(survDat, km, cancerType)
}
dev.off()

# plot forest
survPvaluePlot <- survPvalue
survPvaluePlot$CancerType <- gsub("TCGA_","",survPvaluePlot$CancerType)
survPvaluePlot$Pvalue <- as.numeric(survPvaluePlot$Pvalue)
survPvaluePlot$HR <- as.numeric(survPvaluePlot$HR)
survPvaluePlot$CI_lower <- as.numeric(survPvaluePlot$CI_lower)
survPvaluePlot$CI_upper <- as.numeric(survPvaluePlot$CI_upper)
survPvaluePlot$HR.log <- log2(as.numeric(survPvaluePlot$HR))
survPvaluePlot$CI_lower.log <- log2(as.numeric(survPvaluePlot$CI_lower))
survPvaluePlot$CI_upper.log <- log2(as.numeric(survPvaluePlot$CI_upper))

# cal %
p <- 100*as.numeric(survPvaluePlot$permCAF)/(as.numeric(survPvaluePlot$permCAF)+as.numeric(survPvaluePlot$restCAF))
r <- 100*as.numeric(survPvaluePlot$restCAF)/(as.numeric(survPvaluePlot$permCAF)+as.numeric(survPvaluePlot$restCAF))
survPvaluePlot$permCAF <- paste0(survPvaluePlot$permCAF," (",round(p,0),"%)")
survPvaluePlot$restCAF <- paste0(survPvaluePlot$restCAF," (",round(r,0),"%)")

# derive estimate label
#wrangle results into pre-plotting table form
survPvaluePlot <- survPvaluePlot |>
  # round estimates and 95% CIs to 2 decimal places for journal specifications
  mutate(across(
    c(HR, CI_lower, CI_upper),
    ~ str_pad(
      round(.x, 2),
      width = 4,
      pad = "0",
      side = "right"
    )
  ),
  # add an "-" between HR estimate confidence intervals
  Estimate = paste0(HR, " (", CI_lower, "-", CI_upper, ")")) |>
  # round p-values to 3 decimal places, except in cases where p < .001
  mutate(Pvalue = case_when(
    Pvalue < .001 ~ "<0.001",
    round(Pvalue, 3) == .05 ~ as.character(round(Pvalue,3)),
    Pvalue < .01 ~ str_pad( # if less than .01, go one more decimal place
      as.character(round(Pvalue, 3)),
      width = 5,
      pad = "0",
      side = "right"
    ),
    TRUE ~ str_pad( # otherwise just round to 2 decimal places and pad string so that .2 reads as 0.20
      as.character(round(Pvalue, 3)),
      width = 5,
      pad = "0",
      side = "right"
    )
  )) |>
  # add a row of data that are actually column names which will be shown on the plot in the next step
  bind_rows(
    data.frame(
      CancerType = "Cancer",
      Estimate = "HR (95% CI)",
      permCAF = "permCAF",
      restCAF = "restCAF",
      CI_lower = "",
      CI_upper = "",
      Pvalue = "P-value"
    )
  ) #|>
  #mutate(model = fct_rev(fct_relevel(CancerType, "CancerType")))

#
survPvaluePlot$Estimate[which(survPvaluePlot$CancerType %in% c("TGCT","PRAD"))] <- "N.A."

p_mid <- 
  survPvaluePlot |>
  ggplot(aes(y = reorder(CancerType, HR.log))) + 
  theme_classic() +
  geom_point(aes(x=HR.log), shape=15, size=3) +
  geom_linerange(aes(xmin=CI_lower.log, xmax=CI_upper.log)) +
  geom_vline(xintercept = 0,linetype="dashed") +
  labs(x="Log Hazard Ratio", y="") +
  coord_cartesian(ylim=c(1,34), xlim=c(-4, 6)) +
  #annotate("text", x = -2, y = 35, label = "permCAF better") +
  #annotate("text", x = 2, y = 35, label = "permCAF worse") + 
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank())

p_left <-
  survPvaluePlot  |>
  ggplot(aes(y = reorder(CancerType, HR.log))) +
  geom_text(aes(x = 0, label = CancerType), hjust = 0, fontface = "bold") +
  #geom_text(aes(x = 1, label = Estimate), hjust = 0, fontface = ifelse(survPvaluePlot$Estimate == "HR (95% CI)", "bold", "plain")) +
  geom_text(aes(x = 0.8, label = permCAF), hjust = 0, fontface = ifelse(survPvaluePlot$permCAF == "permCAF", "bold", "plain")) +
  geom_text(aes(x = 1.8, label = restCAF), hjust = 0, fontface = ifelse(survPvaluePlot$restCAF == "restCAF", "bold", "plain")) +
  theme_void() +
  coord_cartesian(ylim=c(1,34), xlim = c(0, 4))

p_right <-
  survPvaluePlot  |>
  ggplot() +
  geom_text(
    aes(x = 0, y = reorder(CancerType, HR.log), label = Pvalue),
    hjust = 0,
    fontface = ifelse(survPvaluePlot$Pvalue == "P-value", "bold", "plain")
  ) +
  theme_void() +
  coord_cartesian(ylim=c(1,34))

# combine
layout <- c(
  area(t = 0, l = 0, b = 34, r = 4), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 0, l = 4, b = 34, r = 5), # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
  area(t = 0, l = 5, b = 34, r = 6) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)
# final plot arrangement
p_left + p_mid + p_right + plot_layout(design = layout)

ggsave("../final_figures/pancan_forest.pdf", width=6, height=6)

