# IMML
An interpretable multimodal machine learning framework that utilizes prior biological knowledge to study complex multifactorial diseases.

## User guide
Currently IMML is implemented as an **R** package. However it will be implemented as a command line tool for greater compatibility. The program's behavior will be customizable trough a config file. For ease of use the following directory structure can be used to avoid specifying the location of each input file.
```
. 
├── input
│   ├── sample_IDs.rds
│   ├── phenotypes.rds
│   ├── ...
│   └── clinical_data.rds
└── output
    ├── clinical_selected.rds
    ├── ...
    ├── final_model.rds
    └── plotting.html
```



### Feature Selection

#### targeted assay
#### untargeted assay
#### clinical data

### Modeltraining
The model incorporates multi-view learning. A critical assumption of multi-view learning, however, is that the single-view models should be independent. This assumption is often violated in complex metabolic diseases, as there is a high level of redundancy and correlations amongst feature layers. Nevertheless, multi-view learning has proven to be superior compared to models leveraging concatenated feature space in crowd-sourced computational challenges. 

#### Authors
Ulrich Asemann & Wilhelm Glaas
