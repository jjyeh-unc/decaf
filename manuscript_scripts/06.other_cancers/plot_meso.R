library(stringr)
library(openxlsx)
library(survival)
library(survminer)
library(survivalAnalysis)
library(ggpubr)

# load
pathCall <- read.xlsx("../data/TCGA_pancan_clinical_caf11523_AI_231114.xlsx")
#clinicalDat <- read.xlsx("../data/TCGA_RNAseq/TCGA-CDR-SupplementalTableS1.xlsx")
names(pathCall)[which(names(pathCall) %in% c("CAPER.score","CAPER","CAPER_graded"))] <- c("DeCAF_prob","DeCAF","DeCAF_graded")
pathCall$DeCAF <- gsub("myCAF","permCAF",pathCall$DeCAF)
pathCall$DeCAF <- gsub("iCAF","restCAF",pathCall$DeCAF)
pathCall$DeCAF_graded <- gsub("myCAF","permCAF",pathCall$DeCAF_graded)
pathCall$DeCAF_graded <- gsub("iCAF","restCAF",pathCall$DeCAF_graded)

pathCall.meso <- pathCall[which(pathCall$type %in% "MESO"),]

pdf("../figure_DeCAF/MESO_barplot_boxplot.pdf")
# barplot
datTmp <- as.data.frame(table(pathCall.meso[,c("DeCAF","histological_type")]))
datTmp$histological_type <- gsub(" mesothelioma", "", datTmp$histological_type)
datTmp$histological_type <- gsub(" - NOS", "", datTmp$histological_type)
datStat <- fisher.test(table(pathCall.meso[,c("DeCAF","histological_type")]))

bar_plots <- ggplot(datTmp, aes(x=factor(datTmp$histological_type, levels=c("Epithelioid", "Biphasic", "Diffuse malignant", "Sarcomatoid")), y=Freq, fill=DeCAF)) +
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
  labs(title = "MESO p = 0.001",
       y = "Proportion of patients",
       x = "Histological type"
  )
bar_plots

# boxplot
datTmp <- data.frame(DeCAF_prob = pathCall.meso$DeCAF_prob,
                     DeCAF =  pathCall.meso$DeCAF,
                     pathCall.di = c("Myxoid","Fibrous")[pathCall.meso$`Collagenized(yes.1,.no.0)`+1],
                     pathCall = c("Myxoid","Mixed","Fibrous")[pathCall.meso$`Myx1/mixed2/Coll3`],
                     stringsAsFactors = FALSE)
datTmp <- datTmp[which(!is.na(datTmp$pathCall.di)),]
datTmp$pathCall.di <- factor(datTmp$pathCall.di, levels = c("Myxoid","Fibrous"))
datTmp$pathCall <- factor(datTmp$pathCall, levels = c("Myxoid","Mixed","Fibrous"))

point_size <- 3
boxplot_lwd <- 0.5
color_mapping_boxplot <- c("CR/PR" = "darkgreen", "SD/PD" = "darkred", "permCAF" = "#FE46A5", "restCAF" = "#008080","Basal-like" = "orange", "Classical" = "blue")


crop2<-ggplot(datTmp, aes(x=pathCall.di, y=DeCAF_prob)) +
  geom_boxplot(outlier.shape = NA, lwd=boxplot_lwd)+
  theme_classic(base_size=15)+
  geom_jitter(aes(color=DeCAF), shape=16, fill="gray", alpha=.75, size=point_size, position=position_jitter(0.2), na.rm=TRUE)+
  scale_color_manual(values = color_mapping_boxplot) +
  stat_compare_means( method = "wilcox.test")+
  labs(x=" ", y="DeCAF Probability",title = "MESO")+
  theme_classic(base_size=15)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
        
  )+
  guides(
    color = guide_legend(override.aes = list(shape = 16,size = 3))
  )

if(FALSE) {
crop2<-ggplot(datTmp, aes(x=pathCall.di, y=DeCAF_prob, color=pathCall.di)) +
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
}

crop2

dev.off()
