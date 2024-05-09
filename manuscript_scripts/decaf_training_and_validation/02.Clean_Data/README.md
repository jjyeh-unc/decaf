## Clean Data Scripts

#### This folder contains the script to clean the data.

# Note: Running this is not necessary to apply the DeCAF classifier to future data. This is included only to demonstrate how the figures for the manuscript were created.


The first part of the script creates functions for cleaning the data for each of the data sets.

The second part of the script loads the DeCAF classifier and the function to apply the DeCAF classifier to data.

The third part of the script applies the cleaning functions to all of the data sets. This is used to ensure that all the data has the same format for the DeCAF classifier evaluations used in the manuscript. 
Lastly, this clean data is saved to the data/clean_data folder. 

Note: Yeh Seq is cleaned and saved separately since it was evaluated separately. 

