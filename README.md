# IMML
An interpretable multimodal machine learning framework that utilizes prior biological knowledge to study complex multifactorial diseases.

## Prerequesites

The conda environment needed for the framework can be installed using the provided `IMML_env.yml` file:
```
conda env create -f IMML_env.yml
```


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

### Data input formats

The `sample_IDs.rds` file needs to be a data frame with samples as rows and columns as the available modalities. Row names are universal sample IDs:

|      | Genomics | Transcriptomics | Proteomics | Metabolomics | Methylomics | Clinical |
|------|----------|-----------------|------------|--------------|-------------|----------|
| 1746 | 2874     | 3576            |       4985 |         5979 | NA          | 7041     |
| 1385 | 2934     | 3327            |         NA |         5234 | 6462        | 7788     |
| 1932 | NA       | 3117            |         NA |         5372 | 6478        | 7468     |

The `phenotypes.rds` file needs to be a data frame with samples as rows and columns as investigated phenotype and covariates:

|      | phenotype | cov1 | cov2 | cov3 |
|------|-----------|------|-----:|-----:|
| 1746 | 0         | 1    |   -5 |  0.4 |
| 1385 | 0         | 1    |    2 | -0.2 |
| 1932 | 1         | 0    |    6 | -0.8 |

### Feature Selection

#### targeted assay
#### untargeted assay
#### clinical data

### Modeltraining
The model incorporates multi-view learning. A critical assumption of multi-view learning, however, is that the single-view models should be independent. This assumption is often violated in complex metabolic diseases, as there is a high level of redundancy and correlations amongst feature layers. Nevertheless, multi-view learning has proven to be superior compared to models leveraging concatenated feature space in crowd-sourced computational challenges. 

## Reference

This framework is based on the analysis presented in the following manuscript in preprint: https://doi.org/10.1101/2024.01.04.574164

Phong BH Nguyen _et al_. The Interpretable Multimodal Machine Learning (IMML) framework reveals pathological signatures of distal sensorimotor polyneuropathy. bioRxiv 2024.01.04.574164.

#### Authors
Ulrich Asemann & Wilhelm Glaas
