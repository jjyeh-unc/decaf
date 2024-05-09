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
pdf("../figure_DeCAF/heatmap_vK.pdf")
# for legend

for (rDataName in c("TCGA_PAAD","CPTAC","Dijk","Moffitt_GEO_array","Grunwald","Hayashi","Linehan","Olive","Puleo_array","PACA_AU_array","PACA_AU_seq")) {
  dataSet <- readRDS(paste("../data_DeCAF/",rDataName,".rds",sep=""))
  cafSubtype <- readRDS(paste("../data_DeCAF/",rDataName,".caf_subtype.rds",sep=""))
  
  cafSubtype$SCISSORS_CAF_K2_top25.unscaled.K2$ColSideColors[which(cafSubtype$SCISSORS_CAF_K2_top25.unscaled.K2$ColSideColors %in% "darkgreen")] <- "violetred1"
  cafSubtype$SCISSORS_CAF_K2_top25.unscaled.K2$ColSideColors[which(cafSubtype$SCISSORS_CAF_K2_top25.unscaled.K2$ColSideColors %in% "coral")] <- "turquoise4"
  
  # renames
  names(cafSubtype)[2] <- "SCISSORS_CAF_K2_top25.unscaled.vK"
  cafSubtype$SCISSORS_CAF_K2_top25.unscaled.vK$ColSideColors[which(cafSubtype$SCISSORS_CAF_K2_top25.unscaled.vK$ColSideColors %in% "darkgreen")] <- "violetred1"
  cafSubtype$SCISSORS_CAF_K2_top25.unscaled.vK$ColSideColors[which(cafSubtype$SCISSORS_CAF_K2_top25.unscaled.vK$ColSideColors %in% c("palegreen3","#7FB17F"))] <- "#FF8BC0" #"#AAD6D8"
  
  cafSubtype$SCISSORS_CAF_K2_top25.unscaled.vK$ColSideColors[which(cafSubtype$SCISSORS_CAF_K2_top25.unscaled.vK$ColSideColors %in% "coral")] <- "turquoise4"
  cafSubtype$SCISSORS_CAF_K2_top25.unscaled.vK$ColSideColors[which(cafSubtype$SCISSORS_CAF_K2_top25.unscaled.vK$ColSideColors %in% c("peachpuff1","#FFBFA7"))] <- "#65B6B9" #"#FFBEDC"
  
  ColSideColors <- data.frame(SCISSORS_vK = cafSubtype$SCISSORS_CAF_K2_top25.unscaled.vK$ColSideColors,
                              SCISSORS = cafSubtype$SCISSORS_CAF_K2_top25.unscaled.K2$ColSideColors,
                              Elyada = cafSubtype$Elyada_CAF.unscaled.K2$ColSideColors,
                              MoffittStroma = cafSubtype$MS.unscaled.K2$ColSideColors,
                              Maurer = cafSubtype$Maurer.unscaled.K2$ColSideColors,
                              stringsAsFactors = FALSE)
  
  #Plot_heatmap_CC_CAF(cafSubtype$SCISSORS_CAF_K2_top25.unscaled.K2, "SCISSORS_CAF_K2_top25", ColSideColors, paste0(rDataName," SCISSORS"))
  Plot_heatmap_CC_CAF(cafSubtype$SCISSORS_CAF_K2_top25.unscaled.vK, "SCISSORS_CAF_K2_top25", ColSideColors, paste0(rDataName," SCISSORS"))
  
  
}

plot.new()
legend(xy.coords(x=0,y=.98),
       legend=c("SCISSORS vK","(training/validation labels)","permCAF","Mixed permCAF","Mixed",
                "Mixed restCAF","restCAF","Absent",
                "SCISSORS K2","myCAF-like","iCAF-like",
                "Elyada","iCAF","myCAF",
                "Moffitt stroma","Activated","Normal",
                "Maurer","ECM-rich","Immune-rich"),
       fill=c("white","white","violetred1","#FF8BC0","grey","#65B6B9","turquoise4","Black",
              "white","violetred1","turquoise4",
              "white","darkgreen","coral",
              "white","brown","skyblue",
              "white","purple3","forestgreen"),
       border=FALSE, bty="n",
       x.intersp = 1 , cex = 1, horiz=FALSE)

dev.off()

