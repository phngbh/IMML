#!/usr/bin/env Rscript

path = funr::get_script_path()
source(file.path(path, "feature_selection.R"))
library(optparse)

get_args <- function(
    # Parse arguments passed to the R script
){
  # Define the command line options
  option_list <- list(
    make_option(c("--data"), type="character", default=NULL, help="path to a dataframe of processed data"),
    make_option(c("--target"), type="character", default=NULL, help="path to a dataframe of target variable and (optional) technical variables"),
    make_option(c("--modals_ids"), type="character", default=NULL, help="path to a dataframe of modality IDs"),
    make_option(c("--target_name"), type="character", default=NULL, help="string of target name (must be a column name in target dataframe)"),
    make_option(c("--pathways"), type="character", default=NULL, help="path to a list of methylomics pathways"),
    make_option(c("--dualmap_file"), type="character", default=NULL, help="path to a dualmap file"),
    make_option(c("--DM_threshold"), type="character", default='2.4e-07', help="significance threshold for differential methylation analysis"),
    make_option(c("--seed"), type="integer", default=993, help="random seed"),
    make_option(c("--n_iterations"), type="integer", default=100, help="number of resamples"),
    make_option(c("--p"), type="double", default=0.8, help="percentage used for resampling"),
    make_option(c("--partitioning"), type="character", default=NULL, help="path to data partitioning file"),
    make_option(c("--out_dir"), type="character", default=NULL, help="path to output dir"))
  
  # Parse the command line options
  opt_parser <- OptionParser(option_list=option_list)
  opt <- parse_args(opt_parser)
  return(opt)
}

opt <- get_args()

partitioning = readRDS(opt$partitioning)
setwd(opt$out_dir)
dir.create('Methylomics')
setwd('Methylomics')
methylomics_create_samples(fs_samples = partitioning$featureSelection$Methylomics,
                           n_resamplings = opt$n_iterations, p = opt$p,
                           seed = opt$seed)

# load target and modality data frames
modals_ids = readRDS(opt$modals_ids)
modals_ids = rownames_to_column(modals_ids, var = 'ID')
targets = readRDS(opt$target)
targets = rownames_to_column(targets, var = 'ID')

# create target dataframe for Methylomics IDs
merged = merge(modals_ids, targets)
methyl_pheno = merged[, c('Methylomics', opt$target_name)]
methyl_pheno = methyl_pheno[methyl_pheno$Methylomics %in% partitioning$featureSelection$Methylomics,]
rownames(methyl_pheno) = NULL
methyl_pheno = column_to_rownames(methyl_pheno, var = 'Methylomics')
saveRDS(methyl_pheno, file = 'methylomics_phenotypes.rds')

command = paste('Rscript', file.path(path, 'processing_methylation.R'),
                getwd(),
                'methylomics_phenotypes.rds',
                file.path(getwd(), 'methylomics_samples'),
                opt$data,
                opt$n_iterations,
                opt$pathways,
                opt$dualmap_file,
                opt$DM_threshold,
                opt$target_name,
                opt$out_dir
                )
system(command)