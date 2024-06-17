# IMML framework
An interpretable multimodal machine learning framework that utilizes prior biological knowledge to study complex multifactorial diseases.

## Prerequesites

The conda environment needed for the framework can be installed using the provided `IMML_env.yml` file:
```
conda env create -f IMML_env.yml
```
This framework uses **nextflow** to efficiently run in parallel on the user's system. Please refer to this link on how to install **nextflow**: https://www.nextflow.io/docs/latest/install.html.

For gene-set analysis of GWAS data this framework uses **MAGMA**. Please refer to this link on how to install **MAGMA**: https://cncr.nl/research/magma/. The installed **MAGMA** executable needs to be added to the IMML conda environment ```bin``` folder. E.g.:
```
mv magma ~/miniconda3/envs/IMML/bin
```



## Usage guide
For an in-depth usage guide please refer to the ```vignette.html```.

First run the feature selection workflow on your command line:
```
nextflow ~/.../IMML/workflows/feature_selection.nf -params-file ~/.../IMML/workflows/FS_config.yml -c ~/.../IMML/nextflow.config
```
And then the model training workflow:
```
nextflow ~/.../IMML/workflows/model_training.nf -params-file ~/.../IMML/workflows/MT_config.yml -c ~/.../IMML/nextflow.config
```
For more information on how to configure **nextflow** using the ```nextflow.config``` file, please refer to this link: https://www.nextflow.io/docs/latest/config.html. In order to use the IMML framework this file needs to include the path to the IMML conda environment:
```
conda.enabled = true
process.conda = '~/miniconda3/envs/IMML/'
```
Input file paths and and parameters are specified in the ```FS_config.yml``` and ```MT_config.yml``` files. See the vignette for more information. 


## Reference

This framework is based on the analysis presented in the following manuscript in preprint: https://doi.org/10.1101/2024.01.04.574164

Phong BH Nguyen _et al_. The Interpretable Multimodal Machine Learning (IMML) framework reveals pathological signatures of distal sensorimotor polyneuropathy. bioRxiv 2024.01.04.574164.

### Authors
Ulrich Asemann, Wilhelm Glaas & Phong BH Nguyen
