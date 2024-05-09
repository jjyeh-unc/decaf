############################## Functions and libraries ##############################
# set dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load libraries
library(ConsensusClusterPlus) # R3.xx and R4.xx have different versions of ConsensusClusterPlus
library("RColorBrewer")
library(openxlsx)
library(stringr)

# load functions
file.sources <- list.files("../R/R/",pattern="*.R")
file.sources <- paste("../R/R/", file.sources, sep="")
sapply(file.sources, source)

# load PurIST
load("../R/PurIST/fitteds_public_2019-02-12.Rdata")
source("../R/PurIST/functions.R")

# load subtype info
### This is a combined subtype object with
### subtypeColList, subtypeGeneList, subtypeList and schemaList
load("../data/cmbSubtypes.RData")
print("Subtype schemas available for use:")
print(schemaList)

############################## Subtyping ##############################
if(FALSE) {
load("../data/TCGA_PAAD_plus.RData")
dataSet <- TCGA_PAAD_plus
rDataName <- "TCGA_PAAD"
# save TME calls
subTME <- read.xlsx("../Grunwald/TCGA_PAAD_stroma-review_FINAL.xlsx")
tmpSplit <- data.frame(str_split_fixed(subTME$sample, "-", 5))
subTME$barcode <- apply(tmpSplit[,c(1:4)], 1, paste, collapse = "-")
subTME <- subTME[match(dataSet$sampInfo$Tumor.Sample.ID, subTME$barcode),]
subTME$predom.stroma.type.per.patient <- factor(subTME$predom.stroma.type.per.patient, levels = c("deserted","inter","reactive","excluded","no slides"))
subTME$`Stromal.amount.(max.10,.10to25,.25to50,.50to80,.min.80)`[which(subTME$`Stromal.amount.(max.10,.10to25,.25to50,.50to80,.min.80)` %in% ">80 ")] <- ">80"
subTME$`Stromal.amount.(max.10,.10to25,.25to50,.50to80,.min.80)` <- factor(subTME$`Stromal.amount.(max.10,.10to25,.25to50,.50to80,.min.80)`, 
                                                                           levels = c("<10", "10to25","25to50","50to80",">80"))
dataSet$subTME <- subTME
# save slide #
dataSet$sampInfo$ImageID.FOR_STATS <- NA
bobby.1 <- read.xlsx("../data/JJY_R5373_20160226slides_Combined_data_BRM_20160509_FOR STATS.xlsx")
tmpSplit <- data.frame(str_split_fixed(bobby.1$`File.Name.(Case)`, "-", 5))
bobby.1$barcode <- apply(tmpSplit[,c(1:4)], 1, paste, collapse = "-")
dataSet$sampInfo$ImageID.FOR_STATS <- bobby.1$`Image.ID.(values)`[match(dataSet$sampInfo$Tumor.Sample.ID, bobby.1$barcode)]

dataSet$sampInfo$ImageID.compileddata3 <- NA
bobby.2 <- read.xlsx("../data/TCGA_anntoations_1312016_compileddata3.xlsx")
tmpSplit <- data.frame(str_split_fixed(bobby.2$`File.Name.(Case)`, "-", 5))
bobby.2$barcode <- apply(tmpSplit[,c(1:4)], 1, paste, collapse = "-")
dataSet$sampInfo$ImageID.compileddata3 <- bobby.2$`Image.ID.(values)`[match(dataSet$sampInfo$Tumor.Sample.ID, bobby.2$barcode)]

# save
saveRDS(dataSet, file = paste("../data_Dec2023/",rDataName,".rds",sep=""))
}
rDataName <- "TCGA_PAAD"
dataSet <- readRDS("../data_DeCAF/TCGA_PAAD.rds")
sampSub <- which(!is.na(dataSet$sampInfo$Grade))
cafSubtype <- list()
cafSubtype$Subtype <- data.frame(sampID = colnames(dataSet$ex), stringsAsFactors = FALSE)
                 
## PurIST
dataSet <- Call_PurIST(dataSet)
cafSubtype$Subtype["PurIST"] <- dataSet$sampInfo$PurIST
cafSubtype$Subtype["PurIST_graded"] <- dataSet$sampInfo$PurIST_graded
cafSubtype$Subtype["PurIST.prob"] <- dataSet$sampInfo$PurIST.prob

if(FALSE) {
## SCISSORS unscaled top10
tmpSchema <- "SCISSORS_CAF_K2_top10"
tmpRowScale <- FALSE
tmpK <- 2
dataSet <- Call_consensus(dataSet, tmpSchema, sampSub, 2, tmpK, tmpRowScale, TRUE, FALSE, "km","euclidean", "R-4.X.X")
tmpSubtypeLab <- c("myCAF","iCAF")
dataSet <- Relabel_subtype(dataSet, tmpSubtypeLab, tmpSchema, tmpK, tmpRowScale)
cafSubtype <- Save_cafSubtype_v2(cafSubtype, dataSet, tmpSchema, tmpK, tmpRowScale, "K2")
}

## SCISSORS unscaled top25
tmpSchema <- "SCISSORS_CAF_K2_top25"
tmpRowScale <- FALSE
tmpK <- 7
dataSet <- Call_consensus(dataSet, tmpSchema, sampSub, 2, tmpK, tmpRowScale, TRUE, FALSE, "km","euclidean", "R-4.X.X")
tmpSubtypeLab <- c("Mixed","myCAF","Mixed.iCAF","myCAF","iCAF","Absent","iCAF")
dataSet <- Relabel_subtype(dataSet, tmpSubtypeLab, tmpSchema, tmpK, tmpRowScale)
cafSubtype <- Save_cafSubtype(cafSubtype, dataSet, tmpSchema, tmpK, tmpRowScale, "vK")

## SCISSORS unscaled top25
tmpSchema <- "SCISSORS_CAF_K2_top25"
tmpRowScale <- FALSE
tmpK <- 2
dataSet <- Call_consensus(dataSet, tmpSchema, sampSub, 2, tmpK, tmpRowScale, TRUE, FALSE, "km","euclidean", "R-4.X.X")
tmpSubtypeLab <- c("myCAF","iCAF")
dataSet <- Relabel_subtype(dataSet, tmpSubtypeLab, tmpSchema, tmpK, tmpRowScale)
cafSubtype <- Save_cafSubtype(cafSubtype, dataSet, tmpSchema, tmpK, tmpRowScale, "K2")

## Elyada unscaled
tmpSchema <- "Elyada_CAF"
tmpRowScale <- FALSE
tmpK <- 2
dataSet <- Call_consensus(dataSet, tmpSchema, sampSub, 2, tmpK, tmpRowScale, TRUE, FALSE, "km","euclidean", "R-4.X.X")
tmpSubtypeLab <- c("myCAF","iCAF")
dataSet <- Relabel_subtype(dataSet, tmpSubtypeLab, tmpSchema, tmpK, tmpRowScale)
cafSubtype <- Save_cafSubtype(cafSubtype, dataSet, tmpSchema, tmpK, tmpRowScale)

## MS unscaled
tmpSchema <- "MS"
tmpRowScale <- FALSE
tmpK <- 2
dataSet <- Call_consensus(dataSet, tmpSchema, sampSub, 2, tmpK, tmpRowScale, TRUE, FALSE, "km","euclidean", "R-4.X.X")
tmpSubtypeLab <- c("Activated","Normal")
dataSet <- Relabel_subtype(dataSet, tmpSubtypeLab, tmpSchema, tmpK, tmpRowScale)
cafSubtype <- Save_cafSubtype(cafSubtype, dataSet, tmpSchema, tmpK, tmpRowScale)

## Maurer
tmpSchema <- "Maurer"
tmpRowScale <- FALSE
tmpK <- 2
dataSet <- Call_consensus(dataSet, tmpSchema, sampSub, 2, tmpK, tmpRowScale, TRUE, FALSE, "km","euclidean", "R-4.X.X")
tmpSubtypeLab <- c("ECM-rich","Immune-rich")
dataSet <- Relabel_subtype(dataSet, tmpSubtypeLab, tmpSchema, tmpK, tmpRowScale)
cafSubtype <- Save_cafSubtype(cafSubtype, dataSet, tmpSchema, tmpK, tmpRowScale)

## Puleo 
tmpSchema <- "Puleo"
dataSet <- Call_centroid(dataSet = dataSet, 
                         tmpSchema = tmpSchema, # this gene set fit the format
                         sampSub = sampSub,
                         lowVarFlt = FALSE, # filter low vairance genes or not; set this to FALSE first. If you get error for low variance of genes, set to TRUE.
                         distance = "euclidean", # see ?ConsensusClusterPlus() for other options
                         Rversion = "R-4.X.X" # or "R-3.X.X" # different R versions use different clusterAlg (km vs kmdist)
)
cafSubtype$Puleo <- dataSet$Puleo
cafSubtype$Subtype$Puleo <- dataSet$sampInfo$Puleo
cafSubtype$Subtype$Puleo.classifier <- FALSE
cafSubtype$Subtype$Puleo.classifier[which(!is.na(cafSubtype$Subtype$Puleo))] <- TRUE

# save
saveRDS(cafSubtype, file = paste("../data_DeCAF/",rDataName,".caf_subtype.rds",sep = ""))

# plot 
#Plot_merged_heatmap.2023(rDataName, dataSet, cafSubtype, sampSub)
