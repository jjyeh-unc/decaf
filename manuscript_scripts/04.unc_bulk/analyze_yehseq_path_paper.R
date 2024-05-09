
library(stringr)
library(openxlsx)
library(survival)
library(survminer)
library(survivalAnalysis)
library(ggpubr)

# load PurIST
load("../R/PurIST/fitteds_public_2019-02-12.Rdata")
source("../R/PurIST/functions.R")

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

# compare calls

dataSet <- readRDS("yehseq_pdac_pi.decaf_freeze.rds")

## parse clincial ---------------------------------------------------------------------
survDat <- data.frame(sampID = colnames(dataSet$ex), 
                      time = dataSet$clinicalDat$`Overall.Survival.(OR.until.death)`,
                      event = dataSet$clinicalDat$Last.Event,
                      neoadj.tx = dataSet$clinicalDat$Neoadj.Tx,
                      neoadj.tx.clean = dataSet$clinicalDat$neo_reg_clean,
                      M = dataSet$clinicalDat$M,
                      DeCAF = dataSet$subtypeCall$DeCAF,
                      DeCAF_prob = dataSet$subtypeCall$DeCAF.prob,
                      PurIST = dataSet$subtypeCall$PurIST,
                      pathCall = dataSet$pathCall$`Myx1/mixed2/Coll3`,
                      pathCall.di = dataSet$pathCall$`Myxoid(yes.1,.no.0)`,
                      whitelist = FALSE,
                      stringsAsFactors = FALSE)
survDat$event[survDat$event == 500 | 
                 survDat$event ==900 |
                 survDat$event ==600 |
                survDat$event ==650 ] <- 0
survDat$event[survDat$event == 2 | 
                 survDat$event ==3] <- 1
survDat$event <- as.numeric(survDat$event)
survDat$time <- as.numeric(survDat$time)
survDat$DeCAF <- factor(survDat$DeCAF , levels = c("permCAF","restCAF"))
survDat$PurIST <- factor(survDat$PurIST, levels = c("Basal-like","Classical"))
survDat$pathCall <- c("Myx","Mixed","Coll")[survDat$pathCall]
survDat$pathCall <- factor(survDat$pathCall, levels = c("Myx","Mixed","Coll"))
survDat$pathCall.di <- c("Coll","Myx")[as.numeric(survDat$pathCall.di)+1]
survDat$pathCall.di <- factor(survDat$pathCall.di, levels = c("Myx","Coll"))
survDat0 <- survDat


# survival

pdf("../figure_DeCAF/yehseq_path.pdf")

# pathCall --------------------------------------------------------------------
survDat <- survDat0
km <- with(survDat, Surv(time,event))
p = coxph(km ~ pathCall, data = survDat)
bic = round(BIC(p),3)
#hr <- round(1/summary(p)$coefficients[2],3)
km_fit <- survfit(km ~ pathCall, data = survDat, type = "kaplan-meier")
p <- ggsurvplot(font.legend = 16, 
                font.title = 16,        
                font.x = 16,  
                font.y = 16, 
                font.tickslab = 16,
                risk.table.fontsize = 7,
                risk.table.height = 0.35,
                surv.median.line = "hv", 
                size = 0.3, 
                censor.size=7, 
                pval.size = 8,
                pval.coord = c(96, 0.75),
                xlim = c(0,144),
                km_fit, conf.int = F, pval = T,
                          legend.title="",break.time.by = 12,
                          legend.labs=c("Myxoid","Mixed","Fibrous"),
                          palette = c("tomato1","khaki2","lightskyblue3"),
                          xlab = "Time (months)", risk.table = T,
                          title = paste("YehSeq pathCall BIC=",bic,sep=""))
print(p)

# pathCall --------------------------------------------------------------------
survDat <- survDat0
km <- with(survDat, Surv(time,event))
p = coxph(km ~ pathCall.di, data = survDat)
bic = round(BIC(p),3)
hr <- round(1/summary(p)$coefficients[2],3)
km_fit <- survfit(km ~ pathCall.di, data = survDat, type = "kaplan-meier")
p <- ggsurvplot(font.legend = 16, 
                font.title = 16,        
                font.x = 16,  
                font.y = 16, 
                font.tickslab = 16,
                risk.table.fontsize = 7,
                risk.table.height = 0.35,
                surv.median.line = "hv", 
                size = 0.3, 
                censor.size=7, 
                pval.size = 8,
                pval.coord = c(96, 0.75),
                xlim = c(0,144),
                km_fit, conf.int = F, pval = T,
                          legend.title="",break.time.by = 12,
                          legend.labs=c("Myxoid","Fibrous"),
                          palette = c("tomato1","lightskyblue3"),
                          xlab = "Time (months)", risk.table = T,
                          title = paste("YehSeq pathCall HR=", hr, ", BIC=",bic,sep=""))
print(p)


# boxplot --------------------------------------------------------------------------------------
survDat <- survDat0[which(!is.na(survDat0$pathCall.di)), ]

crop2<-ggplot(survDat, aes(x=pathCall.di, y=DeCAF_prob, color=pathCall.di)) +
  geom_boxplot(outlier.shape = NA, lwd=1) +
  scale_color_manual(values=c("tomato1","lightskyblue3")) +
  theme_classic(base_size=15)+
  geom_jitter(aes(color=pathCall.di), shape=16, fill="gray", alpha=.75, size=2, position=position_jitter(0.2), na.rm=TRUE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)) +
  stat_compare_means()
  
crop2

crop2<-ggplot(survDat, aes(x=pathCall, y=DeCAF_prob, color=pathCall)) +
  geom_boxplot(outlier.shape = NA, lwd=1) +
  scale_color_manual(values=c("tomato1","khaki2","lightskyblue3")) +
  theme_classic(base_size=15)+
  geom_jitter(aes(color=pathCall), shape=16, fill="gray", alpha=.75, size=2, position=position_jitter(0.2), na.rm=TRUE)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
           axis.text.y = element_text(size = 25),
           axis.title.x = element_text(size = 25),
           axis.title.y = element_text(size = 25),
           legend.text = element_text(size = 15),
           legend.title = element_text(size = 15)) +
  stat_compare_means()
  
crop2


dev.off()