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
    make_option(c("--seed"), type="integer", default=993, help="random seed"),
    make_option(c("--resampling"), type="logical", default=TRUE, help="Whether should do the analysis in many resamples"),
    make_option(c("--n_iterations"), type="integer", default=100, help="number of resamples (if resampling == TRUE)"),
    make_option(c("--p"), type="double", default=0.8, help="percentage used for resampling"),
    make_option(c("--partitioning"), type="character", default=NULL, help="path to data partitioning file"),
    make_option(c("--out_dir"), type="character", default=NULL, help="path to output directory"))
  
  # Parse the command line options
  opt_parser <- OptionParser(option_list=option_list)
  opt <- parse_args(opt_parser)
  return(opt)
}

opt <- get_args()

partitioning = readRDS(opt$partitioning)

data = readRDS(opt$data)
data = data[row.names(data) %in% partitioning$featureSelection$Clinical, ]

# load target and modality data frames
modals_ids = readRDS(opt$modals_ids)
modals_ids = rownames_to_column(modals_ids, var = 'ID')
modals_ids = modals_ids[, c('ID', 'Clinical')]
targets = readRDS(opt$target)
targets = rownames_to_column(targets, var = 'ID')

# create target dataframe for Clinical IDs
merged = merge(modals_ids, targets)
pheno = merged[merged$Clinical %in% partitioning$featureSelection$Clinical,]
rownames(pheno) = NULL
pheno = column_to_rownames(pheno, var = 'Clinical')
pheno = pheno[-1]

result = featureSelection_clinical(
  data = data,
  target = pheno,
  target_name = opt$target_name,
  seed = opt$seed,
  resampling = opt$resampling,
  n_iterations = opt$n_iterations
)

saveRDS(result, file = file.path(opt$out_dir, paste0('selectionRes_Clinical_', opt$target_name, '.rds')))