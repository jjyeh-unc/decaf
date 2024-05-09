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


dataSet <- readRDS("../data_DeCAF/TCGA_PAAD.rds")
dataSet <- Call_PurISS(dataSet, version = "final")
subTME <- dataSet$subTME
subTME$subTME <- NA
subTME$subTME[which(subTME$predom.stroma.type.per.patient %in% c("deserted"))] <- "Deserted"
subTME$subTME[which(subTME$predom.stroma.type.per.patient %in% c("inter","reactive"))] <- "Reactive/intermediate"
subTME$DeCAF <- dataSet$sampInfo$PurISS.final
subTME$DeCAF <- gsub("iCAF","restCAF",subTME$DeCAF)
subTME$DeCAF <- gsub("myCAF","permCAF",subTME$DeCAF)


datTmp <- as.data.frame(table(subTME[,c("DeCAF","subTME")]))



bar_plots <- ggplot(datTmp, aes(x=factor(datTmp$subTME, levels=c("Reactive/intermediate", "Deserted")), y=Freq, fill=DeCAF)) +
  geom_bar(position="fill", stat = "identity")+
  scale_fill_manual(values = alpha(c("permCAF" = "#FE46A5", "restCAF" = "#008080"))) +
  geom_text(aes(label = paste0(Freq)), position = position_fill(vjust = 0.5),
            col="white", size=6)+
  #geom_text(aes(x = 1.5, y = 1.05,
  #              label = paste("Fisher's Exact, p =", signif(datStat$p.value, digits = 3))),
  #          hjust = 0, vjust = -0.5, color = "black", size = 5, family = "Helvetica")+
  theme_classic(base_size=15)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
        
  )+
  labs(title = "TCGA_PAAD",
       y = "Proportion of patients",
       x = "Grunwald subTME"
  )


pdf("../figure_DeCAF/tcga_grunwald_decaf.pdf")
bar_plots

dev.off()
