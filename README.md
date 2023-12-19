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



### DataPartitioning
Data partitioning is performed using the _phenotypes.rds_ and _sample_IDs.rds_ files in the data_partitioning folder. It creates three lists for model training, feature selection training and feature selection testing respectively. The model training partition consists of patients for which all modalities are available. For more details check the documentation below.

<details>
  <summary> Documentation </summary>
  
  #### Description
  
  Builds a list data partitions for the provided sets of IDs.
  
  
  #### Usage
  
  ```r
  DataPartitioning(
    phenotypeIDs,
    dataIDs,
    phenotype,
    partitioning = 0.8,
    numPartitions = 100,
    k = 5,
    iter = 100,
    seed = 123
  )
  ```
  
  
  #### Arguments
  
  Argument      |Description
  ------------- |----------------
  `phenotypeIDs`     |     A data.frame with samples as rows and sample IDs as row names. Columns are phenotypes of interest. If for a sample no information about a phenotype is available, it has to be indicated by NA.
  `dataIDs`     |     A data.frame with samples as rows and the data modalities as columns. Row names represent sample IDs. It holds the data IDs of a sample for each modality. If for a sample there is no data for a modality, it has to be indicated by NA.
  `phenotype`     |     The name of the column in phenotypeIDs, which will be used in the analysis.
  `partitioning`     |     A value, that determines what percentage of samples are used for training during feature selection and model training.
  `numPartitions`     |     The number of partitions that are build for feature selection and model training.
  `k`     |     The amount of folds used for k-fold cross validation during model training.
  `iter`     |     The amount of iterations for k-fold cross validation during model training.
  `seed`     |     The seed used for random number generation. Using the same seed ensures reproducibility.
  
</details>


### FeatureSelection functions
Feature selection is individually performed on each modality. 

#### Clinical
Feature selection for clinical features is performed using elastic net regression. For more details check the documentation below.

<details>
  <summary> Documentation </summary>
  
#### Description

Feature selection for clinical data using elastic net models.


#### Usage

```r
FsClinical(
  trainIDs,
  testIDs,
  dataIDs,
  phenotypeIDs,
  phenotype,
  clinicalData,
  seed = 123
)
```


#### Arguments

Argument      |Description
------------- |----------------
`trainIDs`     |     Set of training IDs from the `DataPartitioning()` function for the feature selection.
`testIDs`     |     Set of testing IDs from the `DataPartitioning()` function for the feature selection.
`dataIDs`     |     A data.frame with samples as rows and the data modalities as columns. Row names represent sample IDs. It holds the data IDs of a sample for each modality. If for a sample there is no data for a modality, it has to be indicated by NA.
`phenotypeIDs`     |     A data.frame with samples as rows and sample IDs as row names. Columns are phenotypes of interest. If for a sample no information about a phenotype is available, it has to be indicated by NA.
`phenotype`     |     The name of the column in phenotypeIDs, which will be used in the analysis.
`clinicalData`     |     A data.frame with samples as rows and sample IDs as row names. Columns are the clinical variables and are named accordingly.
`seed`     |     The seed used for random number generation. Using the same seed ensures reproducibility.

  
</details>

#### Genomics
Feature selection for genomics data is performed by using SNPs for geneset analysis. For more details check the documentation below.

<details>
  <summary> Documentation </summary>
  
#### Description

Feature selection for genomics data using GSEA.


#### Usage

```r
FsGenomics(
  trainIDs,
  testIDs,
  dataIDs,
  phenotypeIDs,
  phenotype,
  binaryFilePrefix,
  geneLocFile,
  geneSetFile,
  outputPrefix,
  seed = 123
)
```


#### Arguments

Argument      |Description
------------- |----------------
`trainIDs`     |     Set of training IDs from the `data_partitioning()` function for the feature selection.
`testIDs`     |     Set of testing IDs from the `data_partitioning()` function for the feature selection.
`dataIDs`     |     A data.frame with samples as rows and the data modalities as columns. It holds the data IDs of a sample for each modality. If for a sample there is no data for a modality, it has to be indicated by NA.
`phenotypeIDs`     |     A data.frame with samples as rows and sample IDs as row names. Columns are phenotypes of interest. If for a sample no information about a phenotype is available, it has to be indicated by NA.
`phenotype`     |     The name of the column in phenotypeIDs, which will be used in the analysis.
`binaryFilePrefix`     |     The file path to the .bim/.bed files used. Only specify the prefix without the file types.
`geneLocFile`     |     The file path to the gene location file to be used.
`geneSetFile`     |     The file path to the gene set file to be used.
`outputPrefix`     |     The prefix of the output file.
`seed`     |     The seed used for random number generation. Using the same seed ensures reproducibility.

  
</details>

#### Metabolomics
Feature selection for metabolomics is performed using DE and subsequent GSE analyses.  For more details check the documentation below.

<details>
  <summary> Documentation </summary>
  
#### Description

Feature selection for metabolomics data using GSEA.


#### Usage

```r
FsMetabolomics(
  trainIDs,
  testIDs,
  dataIDs,
  phenotypeIDs,
  phenotype,
  metabolomicsData,
  geneAnnotation,
  pathwayList,
  gseMinSize = 5,
  resampling = TRUE,
  includeDE = FALSE,
  gseThreshold = 0.1,
  deThreshold = 0.05,
  verbose = FALSE,
  seed = 123
)
```


#### Arguments

Argument      |Description
------------- |----------------
`trainIDs`     |     Set of training IDs from the `data_partitioning()` function for the feature selection.
`testIDs`     |     Set of testing IDs from the `data_partitioning()` function for the feature selection.
`dataIDs`     |     A data.frame with samples as rows and the data modalities as columns. It holds the data IDs of a sample for each modality. If for a sample there is no data for a modality, it has to be indicated by NA.
`phenotypeIDs`     |     A data.frame with samples as rows and sample IDs as row names. Columns are phenotypes of interest. If for a sample no information about a phenotype is available, it has to be indicated by NA.
`phenotype`     |     The name of the column in phenotypeIDs, which will be used in the analysis.
`metabolomicsData`     |     A table holding the data for the metabolomics.
`geneAnnotation`     |     A data.frame of gene IDs and corresponding gene symbols, has at least 3 columns: "ID", "Symbol" and "EntrezID".
`pathwayList`     |     A list containing the pathways for the metabolomics.
`gseMinSize`     |     The minimum amount of genes in the genesets considered during GSE.
`resampling`     |     A logical value, whether resampling should be performed. Resampling goes through the amount of iterations present in `trainIDs` and `testIDs` .
`includeDE`     |     A logical value, whether results from the DE analysis should be included in the final results.
`gseThreshold`     |     The significance threshold for GSE.
`deThreshold`     |     The significance threshold for DE.
`verbose`     |     A logical value, whether verbose information should be printed to the console.
`seed`     |     The seed used for random number generation. Using the same seed ensures reproducibility.
  
</details>

#### Methylomics
work in progress

#### Proteomics
Feature selection for proteomics is performed using DE and subsequent GSE analyses.  For more details check the documentation below.

<details>
  <summary> Documentation </summary>
  
#### Description

The function for the feature selection for proteomics data.


#### Usage

```r
FsProteomics(
  trainIDs,
  testIDs,
  dataIDs,
  phenotypeIDs,
  phenotype,
  proteomicsData,
  geneAnnotation,
  pathwayList,
  gseMinSize = 5,
  resampling = TRUE,
  includeDE = FALSE,
  gseThreshold = 0.1,
  deThreshold = 0.05,
  verbose = FALSE,
  seed = 123
)
```


#### Arguments

Argument      |Description
------------- |----------------
`trainIDs`     |     Set of training IDs from the `data_partitioning()` function for the feature selection.
`testIDs`     |     Set of testing IDs from the `data_partitioning()` function for the feature selection.
`dataIDs`     |     A data.frame with samples as rows and the data modalities as columns. It holds the data IDs of a sample for each modality. If for a sample there is no data for a modality, it has to be indicated by NA.
`phenotypeIDs`     |     A data.frame with samples as rows and sample IDs as row names. Columns are phenotypes of interest. If for a sample no information about a phenotype is available, it has to be indicated by NA.
`phenotype`     |     The name of the column in phenotypeIDs, which will be used in the analysis.
`proteomicsData`     |     A table holding the data for the proteomics.
`geneAnnotation`     |     A data.frame of gene IDs and corresponding gene symbols, has at least 3 columns: "ID", "Symbol" and "EntrezID".
`pathwayList`     |     A list containing the pathways for the metabolomics.
`gseMinSize`     |     The minimum amount of genes in the genesets considered during GSE. Resampling goes through the amount of iterations present in `trainIDs` and `testIDs` .
`resampling`     |     A logical value, whether resampling should be performed.
`includeDE`     |     A logical value, whether results from the DE analysis should be included in the final results.
`gseThreshold`     |     The significance threshold for GSE.
`deThreshold`     |     The significance threshold for DE.
`verbose`     |     A logical value, whether verbose information should be printed to the console.
`seed`     |     The seed used for random number generation. Using the same seed ensures reproducibility.
  
</details>

#### Transcriptomics
Feature selection for proteomics is performed using DE and subsequent GSE analyses.  For more details check the documentation below.

<details>
  <summary> Documentation </summary>
  
#### Description

Feature selection for transcriptomics data using GSEA.


#### Usage

```r
FsTranscriptomics(
  trainIDs,
  testIDs,
  dataIDs,
  phenotypeIDs,
  phenotype,
  transcriptomicsData,
  geneAnnotation,
  pathwayList,
  gseMinSize = 5,
  resampling = TRUE,
  includeDE = FALSE,
  gseThreshold = 0.1,
  deThreshold = 0.05,
  verbose = FALSE,
  seed = 123
)
```


#### Arguments

Argument      |Description
------------- |----------------
`trainIDs`     |     Set of training IDs from the `DataPartitioning()` function for the feature selection.
`testIDs`     |     Set of testing IDs from the `DataPartitioning()` function for the feature selection.
`dataIDs`     |     A data.frame with samples as rows and the data modalities as columns. It holds the data IDs of a sample for each modality. If for a sample there is no data for a modality, it has to be indicated by NA.
`phenotypeIDs`     |     A data.frame with samples as rows and sample IDs as row names. Columns are phenotypes of interest. If for a sample no information about a phenotype is available, it has to be indicated by NA.
`phenotype`     |     The name of the column in phenotypeIDs, which will be used in the analysis.
`transcriptomicsData`     |     A data.frame with genes as rows and gene IDs as row names. Columns are samples with the sample IDs as column names.
`geneAnnotation`     |     A data.frame of gene IDs and corresponding gene symbols, has at least 3 columns: "ID", "Symbol" and "EntrezID".
`pathwayList`     |     A list of pathways used in GSEA.
`gseMinSize`     |     The minimum amount of genes in the genesets considered during GSE.
`resampling`     |     A logical value, whether resampling should be performed. Resampling goes through the amount of iterations present in `trainIDs` and `testIDs` .
`includeDE`     |     A logical value, whether results from the DE analysis should be included in the final results.
`gseThreshold`     |     The significance threshold for GSE.
`deThreshold`     |     The significance threshold for DE.
`verbose`     |     A logical values, whether verbose information should be printed to the console.
`seed`     |     The seed used for random number generation. Using the same seed ensures reproducibility.
  
</details>

### Modeltraining
The model incorporates multi-view learning. A critical assumption of multi-view learning, however, is that the single-view models should be independent. This assumption is often violated in complex metabolic diseases, as there is a high level of redundancy and correlations amongst feature layers. Nevertheless, multi-view learning has proven to be superior compared to models leveraging concatenated feature space in crowd-sourced computational challenges. 

#### Author
Ulrich Asemann & Wilhelm Glaas
