# DeCAF: determination of permissive and restraining cancer-associated fibroblast subtypes

DeCAF is currently implemented in R. This repo contains:

A. function (decaf_functions.R) and data (decaf_classifier.Rdata) needed to run DeCAF. Users only need them downloaded to run DeCAF on a new dataset.

B. an example dataframe (example.rds).

C. all data and scripts used to train/validate DeCAF, and generate figures/tables associated with the manuscript (manuscript_data_scripts).

Disclaimer: DeCAF is patent pending and access/use is for not-for-profit research only.


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
# Interpretate DeCAF result
```{r}
The predictions object above contains a n by 3 dataframe.
The 1st column is the predicted probability of each sample of belonging to the "permCAF" subtype. 
The 2nd column is the subtype call based on a predicted probability cutoff of 0.5. Greater than 0.5 indicates the permCAF subtype, and less than 0.5 indicated the restCAF subtype. 
The 3rd column is a graded DeCAF call, indicating the confidence of the call (Strong, Likely, Lean).
```

# Print DeCAF TSP genes
```{r}
classifier$TSP

      [,1]      [,2]      
 [1,] "IGFL2"   "CHRDL1"  
 [2,] "NOX4"    "OGN"     
 [3,] "VSNL1"   "PI16"    
 [4,] "BICD1"   "ANK2"    
 [5,] "NPR3"    "ABCA8"   
 [6,] "ETV1"    "TGFBR3"  
 [7,] "ITGA11"  "FBLN5"   
 [8,] "CNIH3"   "SCARA5"  
 [9,] "COL11A1" "KIAA1217"
```

# Print DeCAF parameters
```{r}
classifier$fit$beta

                0.1062
(Intercept) -8.3790291
indmat1      1.9399809
indmat2      2.2970861
indmat3      0.7868817
indmat4      1.4660046
indmat5      1.2197586
indmat6      1.6365770
indmat7      1.5563258
indmat8      1.8657202
indmat9      2.3576040
```
