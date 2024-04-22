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
    make_option(c("--fam_file"), type="character", default=NULL, help="path to the .fam file"),
    make_option(c("--file_prefix"), type="character", default=NULL, help="path and prefix of input files"),
    make_option(c("--geneloc_file"), type="character", default=NULL, help="path to gene location file"),
    make_option(c("--geneset_file"), type="character", default=NULL, help="path to gene set file"),
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
dir.create('Genomics')
setwd('Genomics')
genomics_create_samples(fam_file = opt$fam_file ,fs_samples = partitioning$featureSelection$Genomics,
                           n_resamplings = opt$n_iterations, p = opt$p,
                           seed = opt$seed)


# load target and modality data frames
modals_ids = readRDS(opt$modals_ids)
modals_ids = rownames_to_column(modals_ids, var = 'ID')
targets = readRDS(opt$target)
targets = rownames_to_column(targets, var = 'ID')
# create phenotype file
fam = read.csv(opt$fam_file, sep = ' ')
fam = fam[, c(1,2)]
colnames(fam) = c('FID', 'IID')
merged = merge(modals_ids, targets)
gen_pheno = merged[, c('Genomics', opt$target_name)]
gen_pheno = gen_pheno[gen_pheno$Genomics %in% partitioning$featureSelection$Genomics,]
gen_pheno = merge(fam, gen_pheno, by.x = 'IID', by.y = 'Genomics')[c(2, 1, 3)] # reorder after merge
gen_pheno[,3] = gen_pheno[,3] + 1
gen_pheno = gen_pheno[order(gen_pheno[,3]),]
write.table(gen_pheno, file = 'phenotype_file.txt', sep = ' ', row.names = FALSE, col.names = TRUE, quote = FALSE)

command = paste('bash', file.path(dirname(path), 'magma.sh'),
                getwd(),
                opt$file_prefix,
                opt$geneloc_file,
                file.path(getwd(), 'genomics_samples'),
                file.path(getwd(), 'phenotype_file.txt'),
                opt$geneset_file,
                opt$target_name,
                opt$n_iterations
)
system(command)