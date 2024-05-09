library(ggplot2)
library(ggpubr)
library(tidyr)
library(stats)
library(forcats)
library(plyr)
library(stats)
library(reshape2)
library(plyr)
library(dplyr)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(GGally)
library(ggfortify)
library(gridExtra)
library(rstatix)
library(ggcorrplot)
# library(impute)
library(fgsea)
# library(ggVennDiagram)
gc()

setwd("/work/users/c/h/changfei/CAPER_Paper/replot_others/data")

box_plots <- list()
bar_plots <- list()
point_size <- 3
boxplot_lwd <- 0.5
color_mapping_boxplot <- c("CR/PR" = "darkgreen", "SD/PD" = "darkred", "permCAF" = "#FE46A5", "restCAF" = "#008080","Basal-like" = "orange", "Classical" = "blue")
# color_mapping_boxplot2 <- c("CR/PR" = "darkgreen", "SD/PD" = "darkred")
# color_mapping_boxplot <- c("CR/PR" = "green", "SD/PD" = "red", "permCAF" = "#FE46A5", "restCAF" = "#008080")



# TCGA---------------------------------------------

tcga<-read.table("New.TCGA.Survival.txt", sep="\t", header=T, row.names=1)
TCGA<-tcga
TCGA$Consensus<-factor(TCGA$Consensus, levels=c("Ba/Sq", "Stroma-rich", "LumP", "LumU", "LumNS", "NE-like"))

TCGA.Count.Consensus<-TCGA %>%
  subset(Consensus %in% c("Ba/Sq", "Stroma-rich", "LumP", "LumU", "LumNS", "NE-like")) %>%
  dplyr::count(Consensus, DeCAF)

TCGA.Consensus.stats<-fisher.test(TCGA$DeCAF, TCGA$Consensus,workspace = 20000000)

bar_plots[[1]] <- ggplot(TCGA.Count.Consensus, aes(x=factor(TCGA.Count.Consensus$Consensus, 
                                                            levels=c("Ba/Sq", "Stroma-rich", 
                                                                     "LumP", "LumU", "LumNS", "NE-like")), y=n, fill=DeCAF)) +
  geom_bar(position="fill", stat = "identity")+
  scale_fill_manual(values = c("#FE46A5", "#008080")) +
  geom_text(aes(label = paste0(n)), position = position_fill(vjust = 0.5),
            col="black", size=6)+
  geom_text(aes(x = 1.5, y = 1.05,
                label = paste("Fisher's Exact, p =", signif(TCGA.Consensus.stats$p.value, digits = 3))),
            hjust = 0, vjust = -0.5, color = "black", size = 5, family = "Helvetica")+
  theme_classic(base_size=15)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
        
  )+
  labs(title = "TCGA",
       y = "patients",
       x = " "
  )




compare_sub=list(c("Ba/Sq", "Stroma-rich"),
                 c("Ba/Sq", "LumP"),
                 c("Ba/Sq", "LumU"),
                 c("Ba/Sq", "LumNS"), 
                 c("Ba/Sq", "NE-like"))

## box plot Consensus
box_plots[[1]] <- ggplot(TCGA, aes(x=Consensus, y=DeCAF.prob)) +
  geom_boxplot(outlier.shape = NA, lwd=boxplot_lwd)+
  geom_jitter(aes(color=DeCAF), shape=16, fill="gray", alpha=.75, size=point_size, position=position_jitter(0.2), na.rm=TRUE)+
  scale_color_manual(values = c("#FE46A5", "#008080")) +
    stat_compare_means(method="kruskal.test")+
  labs(x=" ", y="DeCAF Probability",title = "TCGA")+
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


# IMvigor------------------------------------
imvig<-read.table("DataOutput/IMvigor.CAF.20231209.txt", sep="\t", header=T, row.names=1)
IMVIG<-imvig %>%
  subset(Consensus %in% c("Ba/Sq", "Stroma-rich", "LumP", "LumU", "LumNS", "NE-like"))
IMVIG$Consensus<-factor(IMVIG$Consensus, levels=c("Ba/Sq", "Stroma-rich", "LumP", "LumU", "LumNS", "NE-like"))

compare_sub=list(c("Ba/Sq", "Stroma-rich"),
                 c("Ba/Sq", "LumP"),
                 c("Ba/Sq", "LumU"),
                 c("Ba/Sq", "LumNS"), 
                 c("Ba/Sq", "NE-like"))

IMVIG.na<-subset(IMVIG, !(is.na(binaryResponse)))
IMVIG.Count.binaryResponse<-IMVIG.na %>%
  dplyr::count(DeCAF, binaryResponse)
IMVIG.Count.Response.stat<-fisher.test(IMVIG.na$binaryResponse,IMVIG.na$DeCAF )
p_tmp <- round(IMVIG.Count.Response.stat$p,digits = 3)


## box plot Consensus
set.seed(123)
box_plots[[2]] <- ggplot(IMVIG, aes(x=Consensus, y=DeCAF_prob)) +
  geom_boxplot(outlier.shape = NA, lwd=boxplot_lwd)+
  theme_classic(base_size=15)+
  geom_jitter(aes(color=DeCAF), shape=16, fill="gray", alpha=.75, size=point_size, position=position_jitter(0.2), na.rm=TRUE)+
  scale_color_manual(values = color_mapping_boxplot) +
  stat_compare_means(comparisons=compare_sub, method = "wilcox")+
  labs(x=" ", y="DeCAF Probability",title = "IMvigor")+
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


## bar plot Consensus
IMVIG.Count.Consensus<-IMVIG %>%
  subset(Consensus %in% c("Ba/Sq", "Stroma-rich", "LumP", "LumU", "LumNS", "NE-like")) %>%
  dplyr::count(Consensus, DeCAF)
IMVIG.Count.Consensus.stat<-fisher.test(IMVIG$Consensus,IMVIG$DeCAF )

bar_plots[[2]] <- ggplot(IMVIG.Count.Consensus, aes(x=factor(IMVIG.Count.Consensus$Consensus, levels=c("Ba/Sq", "Stroma-rich", "LumP", "LumU", "LumNS", "NE-like")), y=n, fill=DeCAF)) +
  geom_bar(position="fill", stat = "identity")+
  scale_fill_manual(values = c("#FE46A5", "#008080")) +
  geom_text(aes(label = paste0(n)), position = position_fill(vjust = 0.5),
            col="black", size=6)+
  geom_text(aes(x = 1.5, y = 1.05,
                label = paste("Fisher's Exact, p =", signif(IMVIG.Count.Consensus.stat$p.value, digits = 3))),
            hjust = 0, vjust = -0.5, color = "black", size = 5, family = "Helvetica")+
  theme_classic(base_size=15)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
        
  )+
  labs(title = "IMvigor",
       y = "patients",
       x = " "
  )


# UROMOL----------------------------------------------------
uromol<-read.table("DataOutput/UROMOL.CAF.20231209.txt", sep="\t", header=T, row.names=1)

## UROMOL Consensus
UROMOL<-uromol %>%
  subset(Consensus %in% c("Ba/Sq", "Stroma-rich", "LumP", "LumU", "LumNS", "NE-like"))
UROMOL$Consensus<-factor(UROMOL$Consensus, levels=c("Ba/Sq", "Stroma-rich", "LumP", "LumU", "LumNS", "NE-like"))
UROMOL$Stage.Grade <- factor(UROMOL$Staging,levels=c("LG-Ta","LG-T1", "HG-Ta", "HG-T1", "HG-CIS")) 
compare_sub=list(c("Ba/Sq", "Stroma-rich"),
                 c("Ba/Sq", "LumP"),
                 c("Ba/Sq", "LumU"),
                 c("Ba/Sq", "LumNS"))
## box plot Consensus
set.seed(123)
box_plots[[3]] <- ggplot(UROMOL, aes(x=Consensus, y=DeCAF.prob)) + 
  geom_boxplot(outlier.shape = NA, lwd=boxplot_lwd)+
  theme_classic(base_size=15)+
  geom_jitter(aes(color=DeCAF), shape=16, fill="gray", alpha=.75, size=point_size, position=position_jitter(0.2), na.rm=TRUE)+
  scale_color_manual(values = c("#FE46A5", "#008080")) +
    stat_compare_means(method = "kruskal.test")+
  labs(x=" ", y="DeCAF Probability",title = "UROMOL")+
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

set.seed(123)
box_plots[[3]]

## bar plot Consensus
UROMOL.Count.Consensus<-UROMOL %>%
  subset(Consensus %in% c("Ba/Sq", "Stroma-rich", "LumP", "LumU", "LumNS")) %>%
  dplyr::count(Consensus, DeCAF)
UROMOL.Count.Consensus.stat<-fisher.test(UROMOL$Consensus,UROMOL$DeCAF )

bar_plots[[3]] <- ggplot(UROMOL.Count.Consensus, aes(x=factor(UROMOL.Count.Consensus$Consensus, 
                                                              levels=c("Ba/Sq", "Stroma-rich", "LumP", "LumU", "LumNS")), 
                                                              y=n, fill=DeCAF)) +
  geom_bar(position="fill", stat = "identity")+
  scale_fill_manual(values = c("#FE46A5", "#008080")) +
  geom_text(aes(label = paste0(n)), position = position_fill(vjust = 0.5),
            col="black", size=6)+
  geom_text(aes(x = 1.5, y = 1.05,
                label = paste("Fisher's Exact, p =", signif(UROMOL.Count.Consensus.stat$p.value, digits = 3))),
            hjust = 0, vjust = -0.5, color = "black", size = 5, family = "Helvetica")+
  theme_classic(base_size=15)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
        
  )+
  labs(title = "UROMOL",
       y = "patients",
       x = " "
  )

## UROMOL Stage.Grade
### box_plots UROMOL Stage.Grade


UROMOL.1<-subset(UROMOL, !(is.na(Stage.Grade)))
box_plots[[4]] <- ggplot(UROMOL.1, aes(x=factor(Stage.Grade, levels=c("LG-Ta","LG-T1", "HG-Ta", "HG-T1")), y=DeCAF.prob)) + 
  geom_boxplot(outlier.shape = NA, lwd=boxplot_lwd)+
  theme_classic(base_size=15)+
  geom_jitter(aes(color=DeCAF), shape=16, fill="gray", alpha=.75, size=point_size, position=position_jitter(0.2), na.rm=TRUE)+
  scale_color_manual(values = c("#FE46A5", "#008080")) +
    stat_compare_means(method="kruskal.test")+
  labs(x=" ", y="DeCAF Probability",title = "UROMOL")+
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
  )+
  scale_y_continuous(
    breaks = c(0.00, 0.5, 1.0),
    labels = c("0.0", "0.5", "1.0"),
    limits = c(0, 1.05)
  ) 

### bar_plots UROMOL Stage.Grade
UROMOL.Count.Stage.Grade <-UROMOL.1$ %>%
  subset(Staging %in% c("LG-Ta","LG-T1", "HG-Ta", "HG-T1")) %>%
  dplyr::count(Staging, DeCAF)
UROMOL.Count.Stage.Grade.stat <- fisher.test(UROMOL.1$Stage.Grade,UROMOL.1$DeCAF )

bar_plots[[4]] <- ggplot(UROMOL.Count.Stage.Grade, aes(x=factor(Staging, levels=c("LG-Ta","LG-T1", "HG-Ta", "HG-T1")), y=n, fill=DeCAF)) +
  geom_bar(position="fill", stat = "identity")+
  scale_fill_manual(values = c("#FE46A5", "#008080")) +
  geom_text(aes(label = paste0(n)), position = position_fill(vjust = 0.5),
            col="black", size=6)+
  geom_text(aes(x = 1.5, y = 1.05,
                label = paste("Fisher's Exact, p =", signif(UROMOL.Count.Stage.Grade.stat$p.value, digits = 3))),
            hjust = 0, vjust = -0.5, color = "black", size = 5, family = "Helvetica")+
  theme_classic(base_size=15)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
        
  )+
  labs(title = "UROMOL",
       y = "patients",
       x = " "
  )
bar_plots[[4]]

## UROMOL UROMOL
### b0x_plots UROMOL UROMOL
box_plots[[5]] <- ggplot(UROMOL, aes(x=factor(UROMOL, levels=c("Class 1", "Class 2a", "Class 2b", "Class 3")), y=DeCAF.prob)) + 
  geom_boxplot(outlier.shape = NA, lwd=boxplot_lwd)+
  theme_classic(base_size=15)+
  geom_jitter(aes(color=DeCAF), shape=16, fill="gray", alpha=.75, size=point_size, position=position_jitter(0.2), na.rm=TRUE)+
  scale_color_manual(values = c("#FE46A5", "#008080")) +
    stat_compare_means(method="kruskal.test")+
  labs(x=" ", y="DeCAF Probability",title = "UROMOL")+
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
  )+
  scale_y_continuous(
    breaks = c(0.00, 0.5, 1.0),
    labels = c("0.0", "0.5", "1.0"),
    limits = c(0, 1.05)
  ) 

### bar_plots UROMOL UROMOL
UROMOL.Count.UROMOL<-UROMOL %>%
  subset(UROMOL %in% c("Class 1", "Class 2a", "Class 2b", "Class 3")) %>%
  dplyr::count(UROMOL, DeCAF)

UROMOL.Count.UROMOL.stat <- fisher.test(UROMOL$UROMOL,UROMOL$DeCAF )
bar_plots[[5]] <- ggplot(UROMOL.Count.UROMOL, aes(x=UROMOL, y=n, fill=DeCAF)) +
  geom_bar(position="fill", stat = "identity")+
  scale_fill_manual(values = c("#FE46A5", "#008080")) +
  geom_text(aes(label = paste0(n)), position = position_fill(vjust = 0.5),
            col="black", size=6)+
  geom_text(aes(x = 1.5, y = 1.05,
                label = paste("Fisher's Exact, p =", signif(UROMOL.Count.UROMOL.stat$p.value, digits = 3))),
            hjust = 0, vjust = -0.5, color = "black", size = 5, family = "Helvetica")+
  theme_classic(base_size=15)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
        
  )+
  labs(title = "UROMOL",
       y = "patients",
       x = " "
  )

pdf("UROMOL.Stage.Grade.DeCAF.Subtype.v2.pdf", width = 7.5,height = 7)
set.seed(123)
box_plots[4]
bar_plots[4]
dev.off()

pdf("Bladder.Plots.DeCAF.Subtype.grouped.kruskal.test.v2.pdf", width = 28 ,height = 12)
set.seed(123)
grid.arrange(grobs = c(box_plots[c(1,2,3,4,5)], bar_plots[c(1,2,3,4,5)]), ncol = 5)
dev.off()

# DeCAF Calls all--------------------------------------
data <- read_excel("/work/users/c/h/changfei/CAPER_Paper/datasets_call/data/Supplementary Table S5 - clinical and molecular data - public 1.xlsx")

data <- as.data.frame(data)
unique(data$Dataset)
data$Dataset <- factor(data$Dataset,levels = c("TCGA_PAAD",
                                               "CPTAC",
                                               "Moffitt_GEO_array",
                                               "PACA_AU_seq",
                                               "PACA_AU_array",
                                               "Linehan",
                                               "Puleo_array",
                                               "Olive",
                                               "Dijk",
                                               "Grunwald",
                                               "Hayashi") )

data <- data[which(data$PDAC_Pi_whitelist == "TRUE"),]
table(data$PDAC_Pi_whitelist)

## box Plot
box_plots[[6]] <- ggplot(data = data, aes(x = PurIST, y = DeCAF_prob)) +
  geom_boxplot(aes(color = PurIST), fill = "transparent", lwd = boxplot_lwd, width =0.7)+
  theme_classic(base_size=15)+
  geom_jitter(aes(color=DeCAF), shape=16, fill="gray", alpha=.75, size=point_size, position=position_jitter(0.2), na.rm=TRUE)+
  scale_color_manual(values = c("#FE46A5", "#008080")) +
  stat_compare_means(size = 5)+
  labs(x=" ", y="DeCAF Probability",title = " ")+
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
  )+
  scale_y_continuous(
    breaks = c(0.00, 0.5, 1.0),
    labels = c("0.0", "0.5", "1.0"),
    limits = c(0, 1.05)
  ) 

## bar Plot
data.Count.PurIST<-data %>%
  subset(PurIST %in% c("Basal-like", "Classical")) %>%
  dplyr::count(PurIST, DeCAF)

data.Count.PurIST.stat <- fisher.test(data$PurIST,data$DeCAF )

bar_plots[[6]] <- ggplot(data.Count.PurIST, aes(x=PurIST, y=n, fill=DeCAF)) +
  geom_bar(position="fill", stat = "identity", width =0.7)+
  scale_fill_manual(values = alpha(color_mapping_boxplot)) +
  geom_text(aes(label = paste0(n)), position = position_fill(vjust = 0.5),
            col="black", size=6)+
  geom_text(aes(x = 1.5, y = 1.05,
                label = paste("Fisher's Exact, p =", signif(data.Count.PurIST.stat$p.value, digits = 3))),
            hjust = 0.5, vjust = -0.5, color = "black", size = 5, family = "Helvetica")+
  theme_classic(base_size=15)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
        
  )+
  labs(title = " ",
       y = "patients",
       x = " "
  )  
  
  
  


# output pdf---------------------
pdf("/work/users/c/h/changfei/CAPER_Paper/replot_others/plot/DeCAFpaper_replot_others_six.pdf",width = 42,height = 14)
set.seed(123)
grid.arrange(grobs = c(box_plots[1:6], bar_plots[1:6]), ncol = 6)
dev.off()








