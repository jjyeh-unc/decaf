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

plots <- list()
point_size <- 3
boxplot_lwd <- 0.5
color_mapping_boxplot <- c("CR/PR" = "darkgreen", "SD/PD" = "darkred", "permCAF" = "#FE46A5", "restCAF" = "#008080")
# color_mapping_boxplot2 <- c("CR/PR" = "darkgreen", "SD/PD" = "darkred")
# color_mapping_boxplot <- c("CR/PR" = "green", "SD/PD" = "red", "permCAF" = "#FE46A5", "restCAF" = "#008080")


# IMmotion------------------------
rcc<-read.table("IMmotion.CAF.20231209.txt", sep="\t", header=T, row.names=1)
rcc.atezo<-subset(rcc, ARM %in% "Atezo")
rcc.atezo<-subset(rcc.atezo, !(is.na(BestResponse.2)))
rcc.atezo$BestResponse.2 <- ifelse(rcc.atezo$BestResponse.2 =="CR_PR","CR/PR","SD/PD")

rcc.atezo.Count.binaryResponse<-rcc.atezo %>%
  dplyr::count(DeCAF, BestResponse.2)
rcc.atezo.Count.binaryResponse.stat<-fisher.test(rcc.atezo$BestResponse.2,rcc.atezo$DeCAF )

## boxplot
set.seed(123)
plots[[1]] <- ggplot(rcc.atezo, aes(x=BestResponse.2, y=DeCAF_prob)) + 
  geom_boxplot(aes(color=BestResponse.2),outlier.shape = NA, lwd=boxplot_lwd, width =0.7)+
  geom_jitter(aes(color=DeCAF), shape=16, alpha=.75, size=point_size, position=position_jitter(0.2), na.rm=TRUE)+
  stat_compare_means(size = 5,method = "t.test")+
  # stat_compare_means(size = 5)+
  scale_color_manual(
    values = color_mapping_boxplot,
    labels = c("CR/PR" = "Box:CR/PR", "SD/PD" = "Box:SD/PD", "permCAF" = "Point:permCAF", "restCAF" = "Point:restCAF"),
    breaks = c("CR/PR", "SD/PD", "permCAF", "restCAF")
    
  ) +
  labs(x=" ", y="DeCAF Probability",title = "IMmotion")+
  theme_classic(base_size=15)+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
        
  ) +
  guides(
    color = guide_legend(title = NULL),
    fill = guide_legend(title = NULL)
  )

set.seed(123)
plots[[1]]

## Bar plot
plots[[2]] <- ggplot(rcc.atezo.Count.binaryResponse, aes(x= DeCAF, y=n, fill=BestResponse.2)) +
  geom_bar(position="fill", stat = "identity", width = 0.7)+
  scale_fill_manual(values = alpha(color_mapping_boxplot, alpha = 0.5)) +
  geom_text(aes(label = paste0(n)), position = position_fill(vjust = 0.5),
            col="black", size=6)+
  geom_text(aes(x = 1.5, y = 1.05,
                label = paste("Fisher's Exact, p =", signif(rcc.atezo.Count.binaryResponse.stat$p.value, digits = 3))),
            hjust = 0.5, vjust = -0.5, color = "black", size = 5, family = "Helvetica")+
  theme_classic(base_size=15)+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
        
  )+
  labs(title = "IMmotion",
       y = "patients",
       x = " ",
       fill = "Best response")
plots[[2]]


# IMvigor------------------------------------
imvig<-read.table("IMvigor.CAF.20231209.txt", sep="\t", header=T, row.names=1)
IMVIG<-imvig %>%
  subset(Consensus %in% c("Ba/Sq", "Stroma-rich", "LumP", "LumU", "LumNS"))
IMVIG$Consensus<-factor(IMVIG$Consensus, levels=c("Ba/Sq", "Stroma-rich", "LumP", "LumU", "LumNS"))

compare_sub=list(c("Ba/Sq", "Stroma-rich"),
                 c("Ba/Sq", "LumP"),
                 c("Ba/Sq", "LumU"),
                 c("Ba/Sq", "LumNS"))

IMVIG.na<-subset(IMVIG, !(is.na(binaryResponse)))
IMVIG.Count.binaryResponse<-IMVIG.na %>%
  dplyr::count(DeCAF, binaryResponse)
IMVIG.Count.Response.stat<-fisher.test(IMVIG.na$binaryResponse,IMVIG.na$DeCAF )
p_tmp <- round(IMVIG.Count.Response.stat$p,digits = 3)


## box plot
set.seed(123)
plots[[3]] <- ggplot(IMVIG.na, aes(x=binaryResponse, y=DeCAF_prob)) + 
  geom_boxplot(aes(color=binaryResponse),outlier.shape = NA, lwd=boxplot_lwd, width =0.7)+
  geom_jitter(aes(color=DeCAF), shape=16, alpha=.75, size=point_size, position=position_jitter(0.2), na.rm=TRUE)+
  stat_compare_means(size = 5,method = "t.test")+
  # stat_compare_means(size = 5)+
  scale_color_manual(
    values = color_mapping_boxplot,
    labels = c("CR/PR" = "Box:CR/PR", "SD/PD" = "Box:SD/PD", "permCAF" = "Point:permCAF", "restCAF" = "Point:restCAF"),
    # breaks = c("green", "red", "#FE46A5", "#008080")
    breaks = c("CR/PR", "SD/PD", "permCAF", "restCAF")
    
  ) +
  labs(x=" ", y="DeCAF Probability",title = "IMvigor")+
  theme_classic(base_size=15)+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
        
  ) +
  guides(
    color = guide_legend(title = NULL),
    fill = guide_legend(title = NULL)
  )

set.seed(123)
plots[[3]]

## bar plot
plots[[4]] <- ggplot(IMVIG.Count.binaryResponse, aes(x= DeCAF, y=n, fill=binaryResponse)) +
  geom_bar(position="fill", stat = "identity", width = 0.7)+
  scale_fill_manual(values = alpha(color_mapping_boxplot, alpha = 0.5)) +
  geom_text(aes(label = paste0(n)), position = position_fill(vjust = 0.5),
            col="black", size=6)+
  geom_text(aes(x = 1.5, y = 1.05,
                label = paste("Fisher's Exact, p =", signif(IMVIG.Count.Response.stat$p.value, digits = 3))),
            hjust = 0.5, vjust = -0.5, color = "black", size = 5, family = "Helvetica")+
  theme_classic(base_size=15)+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
        
  )+
  labs(title = "IMvigor",
       y = "patients",
       x = " ",
       fill = "Best response")
plots[[4]]



pdf("/work/users/c/h/changfei/CAPER_Paper/replot_others/plot/DeCAFpaper_replot_others_BestResponse.pdf",width = 28,height = 7)
set.seed(123)
grid.arrange(grobs = plots[1:4], ncol = 4)
dev.off()








