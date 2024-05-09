## Train Classifier Script

#### This folder contains the script to train the classifier.

# Note: Running this is not necessary to apply the DeCAF classifier to future data. This is included only to demonstrate how the DeCAF classifier was created.

The first part of the script loads all of the public data sets. This is used to obtain a list of all the genes that are common across the public data sets. The goal of this is to increase the likelihood that the genes included in the classifier will be present in future studies.

The second part of the script creates the function used to train the classifier.

The third part of the script cleans the data used to train the classifier. This includes defining the training set, removing samples following the exclusion criteria, and defining the genes determined by SCISSORS to define myCAF and iCAF.

Lastly, the fourth part of the script creates and saves the classifier.

