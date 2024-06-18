# IMML

An interpretable multimodal machine learning framework that utilizes
prior biological knowledge to study complex multifactorial diseases.

## Usage guide

First run the feature selection workflow on your command line:

    nextflow ~/.../IMML/workflows/feature_selection.nf -params-file ~/.../IMML/workflows/FS_config.yml -c ~/.../IMML/nextflow.config

And then the model training workflow:

    nextflow ~/.../IMML/workflows/model_training.nf -params-file ~/.../IMML/workflows/MT_config.yml -c ~/.../IMML/nextflow.config

For more information on how to configure **nextflow** using the
`nextflow.config` file, please refer to this link:
<https://www.nextflow.io/docs/latest/config.html>. In order to use the
IMML framework this file needs to include the path to the IMML conda
environment:

    conda.enabled = true
    process.conda = '~/miniconda3/envs/IMML/'

Input file paths and and parameters are specified in the `FS_config.yml`
and `MT_config.yml` files. As an example this is part of the
`FS_config.yml` file:

    package_path: /your/path/IMML/
    targets: /your/path/targets.rds
    target_name: your_name
    modals_ids: /your/path/modals_ids.rds
    out_dir: /your/path/out/
    seed: 993

Please specify the path to the IMML package in both config files.
Additionally, for the `MT_config.yml` file specify the path to the
`FS_config.yml` file.

### Feature Selection

The feature selection step is customized to be appropriate for each data
type. In general we implement generalized gene set enrichment analysis
(GSEA) for molecular data and regularised linear model for clinical
data. By default, targeted assays such as metabolomics and protein
panels are selected using global test (Goeman et al, 2004), and
untargeted assays such as transcriptomics and proteomics are selected
using the conventional GSEA (Subramanian et al, 2005). Genomics data is
selected with MAGMA (de Leeuw et al, 2015) and methylomics is selected
with methylGSA (Ren et al, 2019). In practice, users can alter this
setting, depending on the type of data and research questions. Details
of the selection algorithm could be found in the paper.

You can specify the modalities to be used in feature selection with the
`feature_selection` list:

    feature_selection:
      - modality: Clinical
        data: /your/path/clinical_data.rds
        resampling: False
        n_iterations: 100
        train_split: 0.8
      - modality: Genomics
        file_prefix: /your/path/genomics_processed
        geneloc_file: /your/path/example.gene.loc
        geneset_file: /your/path/geneset.txt
        n_iterations: 100
        train_split: 0.8
      - modality: Methylomics
        data: /your/path/methylomics_data.rds
        pathways: /your/path/pathways_methyl.rds
        dualmap_file: /your/path/dualmap450kEID.rda
        DM_threshold: 2.4e-07
        n_iterations: 100
        train_split: 0.8
      - modality: targeted
        name: Metabolomics
        data: /your/path/metabolomics_data.rds
        pathways: /your/path/metab_annot.rds
        msea_fdr: 0.2
        resampling: True
        n_iterations: 100
        train_split: 0.8
      - modality: untargeted
        name: Proteomics
        data: /your/path/proteomics_data.rds
        pathways: /your/path/prot_annot.rds
        resampling: True
        n_iterations: 100
        train_split: 0.8
      - modality: untargeted
        name: Transcriptomics
        data: /your/path/transcriptomics_data.rds
        pathways: /your/path/transc_annot.rds
        gsea_fdr: 0.2
        resampling: True
        n_iterations: 100
        train_split: 0.8

Of the `targeted` and `untargeted` modalities you can specify as many as
needed, as long as all of them have distinct names, excluding
`Clinical`, `Genomics` and `Methylomics`. These three modalities can
only be specified once.

### Model Training

We implement the iterative forward feature selection (FFS) to integrate
multiple data modalities. It is based on independent runs of five-fold
cross-validation. In each run, we randomly sample 80% of the dataset to
perform five-fold cross-validation and the performance was tested with
the remaining 20% data. Within the inner loop, a 5-fold cross-validation
selected the best data modality to add next. For each fold of inner CV,
a customizable model is trained using a performance metric of choice
(AUROC, AUPRC or weighted log loss) and tested on the left out samples
of that iteration. The prediction performance of the model is tested by
predicting on the outer test set (20% samples). In each step, the model
adds the next best data modality based on increased performance until
all data modalities are included. Details of the training algorithm
could be found in the paper.

As with feature selection, the model training can be customized using
the `MT_config.yml` file:

    package_path: /your/path/IMML/
    fs_config_file: /your/path/FS_config.yml

    n_iterations: 100
    integration: FFS
    algorithm: glmnet 
    p_metric: AUPRC
    feature: GSEA
    from_clinical: True

    annotation_dir: /your/path/here/

The `integration` options include `FFS` and `ensemble`. `algorithm`
options include `glmnet`, `rf` (random forest) and `svmRadial` from the
`caret` package. `p_metric` options include `wLogLoss` (weighted log
loss), `AUROC` and `AUPRC`. `feature` options include `GSEA`,
`thresholding` (DE?) and `union` (both).

### Data input formats

The `modals_ids.rds` file needs to be a data frame with samples as rows
and columns as the available modalities. Row names are universal sample
IDs:

    head(modals_ids, n=3)

    ##      Genomics Transcriptomics Proteomics Metabolomics Methylomics Clinical
    ## 1746     2874            3576       4985         5979          NA     7041
    ## 1385     2934            3327         NA         5234        6462     7788
    ## 1932       NA            3117         NA         5372        6478     7468

The `targets.rds` file needs to be a data frame with samples as rows and
columns as the target of interest and additional covariates. The target
name needs to be the same as the `target_name` variable in the
`FS_config.yml` file.

    head(targets, n=3)

    ##      phenotype cov1 cov2
    ## 1746         0    1  0.4
    ## 1385         0    1  0.2
    ## 1932         1    0  0.8
