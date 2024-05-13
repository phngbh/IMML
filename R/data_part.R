#!/usr/bin/env Rscript

path = funr::get_script_path()
source(file.path(path, "feature_selection.R"))
library(optparse)

get_args <- function(
    # Parse arguments passed to the R script
){
  # Define the command line options
  option_list <- list(
    make_option(c("--target"), type="character", default=NULL, help="path to a dataframe of target variable and (optional) technical variables"),
    make_option(c("--modals_ids"), type="character", default=NULL, help="path to a dataframe of modality IDs"),
    make_option(c("--target_name"), type="character", default=NULL, help="string of target name (must be a column name in target dataframe)"),
    make_option(c("--seed"), type="integer", default=993, help="random seed"),
    make_option(c("--k"), type="integer", default=5, help="k for k-fold CV"),
    make_option(c("--n_iterations"), type="integer", default=100, help="number of resamples for model training"),
    make_option(c("--percent"), type="double", default=0.8, help="percentage of training data for feature selection"),
    make_option(c("--n_partitions"), type="integer", default=100, help="number of resamples for feature selection"),
    make_option(c("--out_dir"), type="character", default=NULL, help="path to output directory"))
  
  # Parse the command line options
  opt_parser <- OptionParser(option_list=option_list)
  opt <- parse_args(opt_parser)
  return(opt)
}

opt <- get_args()

target = readRDS(opt$target)
modals_ids = readRDS(opt$modals_ids)

part = data_partitioning(phenotypeIDs = target, dataIDs = modals_ids,
                         phenotype = opt$target_name,
                         partitioning = opt$percent,
                         numPartitions = opt$n_paritions,
                         k = opt$k,
                         iter = opt$n_iterations,
                         seed = opt$seed)

saveRDS(part, file = file.path(opt$out_dir, paste0('data_partition_', opt$target_name, '.rds')))