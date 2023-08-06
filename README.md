# IML
An interpretable multimodal machine learning framework that utilises prior biological knowledge to study complex multifactorial diseases

## User guide
A guide of the current state of the package for continued development.

### DataPartitioning function
For the first function two main files have to be provided. The files are _phenotypes.rds_ and _sample_IDs.rds_. The other parameters can be set as wished, but the default values are set as the norm that is normaly used. More information is provided in the function itself.

### FeatureSelection functions
Feature selection has a general function that is empty. It was planned to be used as a super function that is called with the respective modality instead of each function single handed. This can be implemented in the future.
The results of the functions are currently saved in the environment of the package, so it is possible to check if the results are right.

#### Clinical
For clinical feature selection the training and testing set from the data partitioning is needed to calculate the data. Also a file with the clinical data is needed to calculate the data. Other files needed are the two main files from the data partitioning. They are needed for all of the other functions from the feature selection.

#### Genomics
For genomics the data is loaded directly inside of the function at the moment. The file path needs to be added inside of the function, where the files are on the computer. Also two additional files are used by the function. Those can be added to the load function to make the load process faster.

#### Metabolomics
The metabolomics function needs the files from the data partitioning, the main ones from the beginning and some additional files. More details are given inside the function itself.

#### Methylomics
The methylomics function is close to the metabolomics function in the case of the used files. Look into the function for more details.

#### Proteomics
The proteomics function is also very close to the the metabolomics, with some differences in the additional files. More information is given in the function itself.

#### Transcriptomics
More details are given inside the function itself.

### Modeltraining functions
For the step of the modeltraining, some functions are already written. 
Following functions were already added:
_MakeOverlapSet_
_SelectModeltrainingIDs_
_Standardise_
_BedToDataframeCopyNumber_

The functions are found in the ml3.R file beginning at the code line 1044.
The functions coming after that need to be added in the future for the modeltraining step.
Also more descriptions need to be added to these functions, when their functionlity is understood.

### Load function
The load function can be used to load all the needed data sheets into the environment of RStudio. In the function the file paths need to be added, where the computer can find them. The important data will always be loaded. The data for the respective modality will only be loaded, if the parameter of the modality is set to TRUE when the load function is called. This ensures not to load data that is not needed for the current development of the package.

### Other functions
Functions that are used additional inside of other functions:
_CaretLogLoss_
_LogLoss_

These functions are used inside of some of the feature selection functions for additional calculations.
More descriptions have to be added in the future.

### Data folder
The data folder currently holds two files that are loaded by functions at run time. Some calculations of feature selection functions can still be put into such files, so the run time can be reduced and the functionality of the package can be made faster.

### Data
Data that is used to build the package and test the functions need to be provided by Ba Hung Phong Nguyen. 

### Future Steps
To continue with the package, the current step evolves around the modeltraining functions. These need to be implemented in the next step of the project. 
Other goals include adding multi core calculations to increase the speed of the package and reduce the time needed. 
Another goal is to generalize the package at some points more, so it is available for other disease and the possibility for other people to use it.

#### Author
Ulrich Asemann
