## Generate Figures

#### This folder contains the script to recreate the figures from this paper.

# Note: Running this is not necessary to apply the DeCAF classifier to future data. This is included only to demonstrate how the figures for the manuscript were created.


The first part of the scripts create functions for generating the figures.

The second part of the script loads the DeCAF classifier and the clean data from the data/clean_data folder generated using the scripts in the 02.Clean_Data folder. 

The third part of the script applies the data generating functions. The figures are saved to data/results folder. 

The first script generates all of the survival figures. The second script generates figures regarding the accuracy of the DeCAF calls compared to the SCISSORS calls.