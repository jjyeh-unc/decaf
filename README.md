# DeCAF: determination of permissive and restraining cancer-associated fibroblast subtypes

- DeCAF is currently implemented in R.
- R version 4.3.3 was used for development/testing.
- R version 4.0.0 or newer is highly recommended.
- Operating systems that are compatible with R version 4.0.0 or newer are required.
- Non-standard hardware is not required.
- Installation time is minimal by donwloading just decaf_functions.R, decaf_classifier.Rdata and example.rds.
- A typical runtime for a dataset with around 150 samples is less than 5 seconds.

This repo contains:

A. function (decaf_functions.R) and data (decaf_classifier.Rdata) needed to run DeCAF, as well as an example dataframe (example.rds).

B. newly generated bulk RNAseq (N=129, UNC-bulk) and scRNAseq (N=6, UNC-sc) data.

C. scripts and other parsed public data used to train/validate DeCAF, and generate figures/tables associated with the manuscript (manuscript_data_scripts).

Disclaimer: DeCAF is patent pending and access/use is for not-for-profit research only.


# Installation
Download decaf_functions.R, decaf_classifier.Rdata and example.rds in the working directory.

# Run DeCAF
```{r}
# load data objects and source functions
load("decaf_classifier.Rdata")
source("decaf_functions.R")
dat = readRDS("example.rds")

# Extract classifier 
classifier = decaf_classifier$classifier2

# apply classifier 
predictions = apply_decaf(data = dat, classifier = classifier)
```

# Interpretate DeCAF output
```{r}
# The predictions object above contains a n (samples) by 3 dataframe.
# The 1st column is the predicted probability of each sample of belonging to the "permCAF" subtype. 
# The 2nd column is the subtype call based on a predicted probability cutoff of 0.5. Greater than 0.5 indicates the permCAF subtype, and less than 0.5 indicated the restCAF subtype. 
# The 3rd column is a graded DeCAF call, indicating the confidence of the call (Strong, Likely, Lean).

print(predictions)

#                  DeCAF_prob   DeCAF   DeCAF_graded
# TCGA.2L.AAQE.01A 0.51673518 permCAF   Lean permCAF
# TCGA.XD.AAUL.01A 0.99882687 permCAF Strong permCAF
# TCGA.3A.A9IN.01A 0.00117836 restCAF Strong restCAF
# TCGA.3A.A9IS.01A 0.44116428 restCAF   Lean restCAF
# TCGA.2L.AAQJ.01A 0.01071204 restCAF Strong restCAF
# TCGA.2L.AAQI.01A 0.06761943 restCAF Strong restCAF
# TCGA.3A.A9IB.01A 0.70824126 permCAF Likely permCAF
# TCGA.3A.A9IU.01A 0.99603845 permCAF Strong permCAF
# TCGA.2L.AAQM.01A 0.51673518 permCAF   Lean permCAF
# TCGA.FB.AAPS.01A 0.98306257 permCAF Strong permCAF
```

# Print DeCAF TSP genes
```{r}
classifier$TSP

 #      [,1]      [,2]      
 # [1,] "IGFL2"   "CHRDL1"  
 # [2,] "NOX4"    "OGN"     
 # [3,] "VSNL1"   "PI16"    
 # [4,] "BICD1"   "ANK2"    
 # [5,] "NPR3"    "ABCA8"   
 # [6,] "ETV1"    "TGFBR3"  
 # [7,] "ITGA11"  "FBLN5"   
 # [8,] "CNIH3"   "SCARA5"  
 # [9,] "COL11A1" "KIAA1217"
```

# Print DeCAF parameters
```{r}
classifier$fit$beta

#                0.1062
# (Intercept) -8.3790291
# indmat1      1.9399809
# indmat2      2.2970861
# indmat3      0.7868817
# indmat4      1.4660046
# indmat5      1.2197586
# indmat6      1.6365770
# indmat7      1.5563258
# indmat8      1.8657202
# indmat9      2.3576040
```


# DeCAF functionality description (pseudocode)
```{r}

# This is the main function to implement the DeCAF classifier.
# This function is associated with Fig. 1a in the manuscript.

apply_decaf = function(data, classifier){ 
  
  # "data": a dataframe with unique official gene symbols as rownames
  # "classifier": the DeCAF classifier containing essential objects
  
  ## Extract Gene Games
  genes = rownames(data)
  
  ## Extract Classifier 
  fit = classifier$fit
  if(is.null(fit$beta)) "Classifier Does Not Have Coefficients Assigned to beta"
  
  ## Keep only gene info for genes in classifier 
  data1 =  data[genes %in% classifier$TSPs,]
  rnames = rownames(data1)
  data1 = matrix(as.numeric(as.matrix(data1)), ncol = ncol(data1))
  rownames(data1) = rnames
  if(nrow(data1) != length(unique(classifier$TSPs))){
    print(classifier$TSPs[!classifier$TSPs %in% genes])
    stop("genes missing")
  }
  
  ## See which of gene pair has higher expression
  indmat = matrix(-1, ncol(data1), nrow(classifier$TSPs))
  for(i in 1:nrow(classifier$TSPs)){
    p1 = which(rownames(data1) == classifier$TSPs[i,1])
    p2 = which(rownames(data1) == classifier$TSPs[i,2])
    indmat[,i] = (data1[p1,] > data1[p2,])^2
  }
  
  ## Calculate probability of permCAF
  X=cbind(rep(1, nrow(indmat)), indmat)
  trainingPrediction = exp(X%*%c(fit$beta))/(1+exp(X%*%c(fit$beta)))
  
  ## Obtain DeCAF subtype
  classification = c("restCAF","permCAF")[(trainingPrediction >= 0.5)^2 + 1]
  
  ## Obtain Grade of DeCAF Subtype
  guess = rep(1, length(trainingPrediction))
  guess[trainingPrediction < .1] = "Strong restCAF"
  guess[trainingPrediction >= .1 & trainingPrediction < .4] = "Likely restCAF"
  guess[trainingPrediction >= .4 & trainingPrediction < .5] = "Lean restCAF"
  guess[trainingPrediction >= .5 & trainingPrediction < .6] = "Lean permCAF"
  guess[trainingPrediction >= .6 & trainingPrediction < .9] = "Likely permCAF"
  guess[trainingPrediction >= .9 ] = "Strong permCAF"
  
  ## Put results together into dataframe
  final = data.frame(DeCAF_prob = trainingPrediction, DeCAF = classification, 
                     DeCAF_graded = guess)
  rownames(final) = make.names(colnames(data), unique = any(table(colnames(data)) > 1) )
  
  return(final)
}
```
