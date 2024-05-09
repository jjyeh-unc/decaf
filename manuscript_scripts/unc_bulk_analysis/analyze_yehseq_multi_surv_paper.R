
library(stringr)
library(openxlsx)
library(survival)
library(survminer)
library(survivalAnalysis)
#library(forestplot)
library(magrittr)
library(dplyr)

dataSet <- readRDS("yehseq_pdac_pi.decaf_freeze.rds")
clinicalDatParsed <- readRDS("yehseq_pdac_pi.decaf_freeze_clinicalParsed_111723.rds")
clinicalDatParsed <- clinicalDatParsed$clinicalDatParsed

## parse clincial ---------------------------------------------------------------------
survDat <- data.frame(sampID = colnames(dataSet$ex), 
                      time = dataSet$clinicalDat$`Overall.Survival.(OR.until.death)`,
                      event = dataSet$clinicalDat$Last.Event,
                      #Neoadj.tx = dataSet$clinicalDat$Neoadj.Tx,
                      #neoadj.tx.clean = dataSet$clinicalDat$neo_reg_clean,
                      M = dataSet$clinicalDat$M,
                      Differentiation = clinicalDatParsed$differentiation_clean_downstaged,
                      T = clinicalDatParsed$T.impute.from.TumorSize.8thEd,
                      N = clinicalDatParsed$N.impute.from.PosLN.8thEd,
                      Lymphovascular.Invasion = clinicalDatParsed$Lymphovascular.Invasion.clean,
                      Margin = as.character(clinicalDatParsed$`Margin.(close.<1mm)`),
                      Stage = clinicalDatParsed$Stage.impute.TN.8thEd,
                      #Adj.tx = clinicalDatParsed$Adj_Tx.clean,
                      #preop_ca199 = clinicalDatParsed$preop_ca199,
                      #pathCall = c("Myx","Mixed","Col")[dataSet$pathCall$`Myx1/mixed2/Coll3`],
                      #pathCall.di = c("Col","Myx")[dataSet$pathCall$`Myxoid(yes.1,.no.0)`+1],
                      DeCAF = dataSet$subtypeCall$DeCAF,
                      PurIST = dataSet$subtypeCall$PurIST,
                      whitelist = FALSE,
                      stringsAsFactors = FALSE)
#survDat$Neoadj.tx[which(survDat$Neoadj.tx %in% "no")] <- "No"
survDat$event[survDat$event == 500 | 
                 survDat$event ==900 |
                 survDat$event ==600 |
                survDat$event ==650 ] <- 0
survDat$event[survDat$event == 2 | 
                 survDat$event ==3] <- 1
survDat$event <- as.numeric(survDat$event)
survDat$time <- as.numeric(survDat$time)
survDat$DeCAF <- factor(survDat$DeCAF, levels = c("restCAF","permCAF"))
survDat$PurIST <- factor(survDat$PurIST, levels = c("Classical","Basal-like"))
survDat$Differentiation <- factor(survDat$Differentiation, levels = c("Well","Moderate","Poor"))
survDat$Stage <- factor(survDat$Stage, levels = c("I","II","III"))
survDat$Margin[which(survDat$Margin %in% c("Close","Negative"))] <- "Close/Negative"
survDat0 <- survDat

pdf("yehsesq_multi_surv.pdf")
# +differentiation
survDat <- survDat0[complete.cases(survDat0), ]

fit <- coxph(formula = Surv(time, event) ~ DeCAF+PurIST+Stage+Differentiation+Margin+Lymphovascular.Invasion, data = survDat)
p <- ggforest(fit, 
              main = "+Differentiation",
              fontsize = 1,
              cpositions = c(0.02, 0.13, 0.32)) 
print(p)

# -differentiation
survDat <- survDat0[, -which(names(survDat0) %in% "Differentiation")]
survDat <- survDat[complete.cases(survDat), ]

fit <- coxph(formula = Surv(time, event) ~ DeCAF+PurIST+Stage+Margin+Lymphovascular.Invasion, data = survDat)
p <- ggforest(fit, 
              main = "-Differentiation",
              fontsize = 1,
              cpositions = c(0.02, 0.13, 0.32)) 
print(p)

# Close/Negative
survDat <- survDat0[survDat0$Margin == "Close/Negative",]
survDat <- survDat[complete.cases(survDat), ]

fit <- coxph(formula = Surv(time, event) ~ DeCAF+PurIST+Stage+Differentiation+Lymphovascular.Invasion, data = survDat)
p <- ggforest(fit, 
              main = "Margin = close or negative, +Differentiation",
              fontsize = 1,
              cpositions = c(0.02, 0.13, 0.32)) 
print(p)

# Close/Negative
survDat <- survDat0[survDat0$Margin == "Close/Negative",]
survDat <- survDat[complete.cases(survDat), ]

fit <- coxph(formula = Surv(time, event) ~ DeCAF+PurIST+Stage+Lymphovascular.Invasion, data = survDat)
p <- ggforest(fit, 
              main = "Margin = close or negative, -Differentiation",
              fontsize = 1,
              cpositions = c(0.02, 0.13, 0.32)) 
print(p)

dev.off()


