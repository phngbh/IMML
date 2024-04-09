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
The feature selection step is customized to be appropriate for each data type. In general we implement generalized gene set enrichment analysis (GSEA) for molecular data and regularised linear model for clinical data. By default, targeted assays such as metabolomics and protein panels are selected using global test (Goeman et al, 2004), and untargeted assays such as transcriptomics and proteomics are selected using the conventional GSEA (Subramanian et al, 2005). Genomics data is selected with MAGMA (de Leeuw et al, 2015) and methylomics is selected with methylGSA (Ren et al, 2019). In practice, users can alter this setting, depending on the type of data and research questions. Details of the selection algorithm could be found in the paper.    

#### targeted assay
#### untargeted assay
#### clinical data

### Modeltraining
We implement the iterative forward feature selection (FFS) to integrate multiple data modalities. It is based on 100 independent runs of five-fold cross-validation. In each run, we randomly sample 80% of the dataset to perform five-fold cross-validation and the performance was tested with the remaining 20% data. Within the inner loop, a 5-fold cross-validation selected the best data modality to add next. For each fold of inner CV, a customizable model is trained using a performance metric of choice (AUROC, AUPRC or weighted log loss) and tested on the left out samples of that iteration. The prediction performance of the model is tested by predicting on the outer test set (20% samples). In each step, the model adds the next best data modality based on increased performance until all data modalities are included. Details of the training algorithm could be found in the paper.

## Reference

This framework is based on the analysis presented in the following manuscript in preprint: https://doi.org/10.1101/2024.01.04.574164

Phong BH Nguyen _et al_. The Interpretable Multimodal Machine Learning (IMML) framework reveals pathological signatures of distal sensorimotor polyneuropathy. bioRxiv 2024.01.04.574164.

#### Authors
Ulrich Asemann, Wilhelm Glaas & Phong BH Nguyen
