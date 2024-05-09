############################## Functions and libraries ##############################
# set dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load libraries
library(ConsensusClusterPlus) # R3.xx and R4.xx have different versions of ConsensusClusterPlus
library("RColorBrewer")

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

# change scissors color
tmpGeneInfo <- subtypeGeneList[[18]]
tmpGeneInfo$Color[which(tmpGeneInfo$Color %in% "darkgreen")] <- "violetred1"
tmpGeneInfo$Color[which(tmpGeneInfo$Color %in% "coral")] <- "turquoise4"
subtypeGeneList[[18]] <- tmpGeneInfo  

############################## plot ##############################
pdf("../figure_DeCAF/heatmap_K2.pdf")
# for legend
plot.new()
legend(xy.coords(x=0,y=.98),
       legend=c("SCISSORS","myCAF-like","iCAF-like",
                "Elyada","iCAF","myCAF",
                "Moffitt stroma","Activated","Normal",
                "Maurer","ECM-rich","Immune-rich"),
       fill=c("white","violetred1","turquoise4",
              "white","darkgreen","coral",
              "white","brown","skyblue",
              "white","purple3","forestgreen"),
       border=FALSE, bty="n",
       x.intersp = 0.5 , cex = 0.35, horiz=TRUE)
for (rDataName in c("TCGA_PAAD","CPTAC","Dijk","Moffitt_GEO_array","Grunwald","Hayashi","Linehan","Olive","Puleo_array","PACA_AU_array","PACA_AU_seq")) {
  dataSet <- readRDS(paste("../data_DeCAF/",rDataName,".rds",sep=""))
  cafSubtype <- readRDS(paste("../data_DeCAF/",rDataName,".caf_subtype.rds",sep=""))
  
  cafSubtype$SCISSORS_CAF_K2_top25.unscaled.K2$ColSideColors[which(cafSubtype$SCISSORS_CAF_K2_top25.unscaled.K2$ColSideColors %in% "darkgreen")] <- "violetred1"
  cafSubtype$SCISSORS_CAF_K2_top25.unscaled.K2$ColSideColors[which(cafSubtype$SCISSORS_CAF_K2_top25.unscaled.K2$ColSideColors %in% "coral")] <- "turquoise4"
  ColSideColors <- data.frame(SCISSORS = cafSubtype$SCISSORS_CAF_K2_top25.unscaled.K2$ColSideColors,
                              Elyada = cafSubtype$Elyada_CAF.unscaled.K2$ColSideColors,
                              MoffittStroma = cafSubtype$MS.unscaled.K2$ColSideColors,
                              Maurer = cafSubtype$Maurer.unscaled.K2$ColSideColors,
                              stringsAsFactors = FALSE)
  
  Plot_heatmap_CC_CAF(cafSubtype$SCISSORS_CAF_K2_top25.unscaled.K2, "SCISSORS_CAF_K2_top25", ColSideColors, paste0(rDataName," SCISSORS"))
  Plot_heatmap_CC_CAF(cafSubtype$Elyada_CAF.unscaled, "Elyada_CAF", ColSideColors, paste0(rDataName," Elyada"))
  Plot_heatmap_CC_CAF(cafSubtype$MS.unscaled.K2, "MS", ColSideColors, paste0(rDataName," MS"))
  Plot_heatmap_CC_CAF(cafSubtype$Maurer.unscaled.K2, "Maurer", ColSideColors, paste0(rDataName," Maurer"))
}

dev.off()

