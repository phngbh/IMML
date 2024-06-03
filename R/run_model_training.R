#!/usr/bin/env Rscript

# libs <- c("Rlib", .libPaths())
# .libPaths(libs)

suppressMessages(library(funr))
suppressMessages(library(yaml))

cat("Load functions and environment\n")
path = funr::get_script_path()
source(file.path(path, "model_training.R"))

cat("Get the argument settings\n")
opt <- get_args()
iter <- opt$iter
integration <- opt$integration
algorithm <- opt$algorithm
p_metric <- opt$p_metric
features <- opt$feature
config_path <- opt$config

config = read_yaml(config_path)
fs_config = config$feature_selection
wdir = config$out_dir

setwd(wdir)
if(!dir.exists("model_results")){ dir.create("model_results") }

cat("Load preprocessed data and necessary files\n")
targets <- readRDS(config$targets)
targets = rownames_to_column(targets, var = 'SampleID')
modals_ids <- readRDS(config$modals_ids)
modals_ids = rownames_to_column(modals_ids, var = 'SampleID')
outcome_name <- config$target_name
preprocessed_data <- list()

for (i in 1:length(fs_config)) {
  modality = fs_config[[i]]
  if (modality$modality == 'Clinical'){
    preprocessed_data[[modality$modality]] = readRDS(modality$data)
  }
  else if (modality$modality == 'Genomics') {
    next
  }
  else if (modality$modality %in% c('untargeted', 'targeted')){
    preprocessed_data[[modality$name]] = t(readRDS(modality$data))
  }
  else if (modality$modality == 'Methylomics'){
    preprocessed_data[[modality$modality]] = t(readRDS(modality$data))
  }
  else {
    stop("Error: unknown modality!")
  }
}

# for (i in names(tmp_dict)){
#   if (file.exists(paste0("processed_",tmp_dict[[i]],"_",outcome_name,".rds"))){
#     preprocessed_data[[i]] <- readRDS(paste0("processed_",tmp_dict[[i]],"_",outcome_name,".rds"))
#     if (!i %in% c("Olink", "Clinical")){
#       preprocessed_data[[i]] <- t(preprocessed_data[[i]])
#     }
#   }
# }

cat("...There are ",length(preprocessed_data), " modalities loaded\n")

cat("Load feature selection results\n")
selection_result <- list()
mod_names = names(preprocessed_data)
for (name in mod_names){
  if (file.exists(paste0("selectionRes_",name,"_",outcome_name,".rds"))){
    selection_result[[name]] <- readRDS(paste0("selectionRes_",name,"_",outcome_name,".rds"))
  }
}
cat("...There are ",length(selection_result), " modality loaded\n")

if (any(!names(preprocessed_data) %in% names(selection_result)) | any(!names(selection_result) %in% names(preprocessed_data))){
  stop("The modality in the provided data and feature selection results do not match")
} else {
  preprocessed_data <- preprocessed_data[names(selection_result)]
}

cat("Make selected data list\n")
selected_data <- make_selected_list(data_list = preprocessed_data, feature_selection_result = selection_result, feature_type = features)
cat("...Free up some memory\n")
rm(preprocessed_data)
rm(selection_result)
gc()
cat("...Add genomics data to selected data list\n")
if (file.exists(paste0("selectionRes_Genomics_", outcome_name,".bed"))){
  selected_data$Genomics <- bed_to_df(paste0("selectionRes_Genomics_", outcome_name,".bed"))
} else {
  cat("......There is no selected genomics data\n")
}

cat("Make input data list\n")
outcome <- targets[,outcome_name]
names(outcome) <- targets$SampleID
outcome <- na.omit(outcome)
cat("There are ", length(outcome), " samples with outcome information in total\n")
input_data <- make_input_list(data_list = selected_data, outcome = outcome, id_table = modals_ids)
# check if shouldn't be overwritten
saveRDS(input_data,file.path("model_results", paste0(outcome_name,"_input_data.rds")))

cat("Make cross validation list\n")
outcome <- outcome[rownames(input_data[[1]])]
cv_list <- make_cv_list(outcome = as.factor(outcome))

cat("Train and evaluate models\n")
outcome <- ifelse(outcome == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
if (opt$from_clinical){
  if (integration == "FFS"){
    res <- fit_forwardSelectFromClinical(data_list = input_data, y = outcome, cv_list = cv_list, p_metric = p_metric, algorithm = algorithm, n = iter)
    saveRDS(res,file.path("model_results", paste0(paste(outcome_name,integration,algorithm,p_metric,features,"fromClinical",iter,sep = "_"),".rds")))
  } else if (integration == "ensemble") {
    res <- fit_ensemble(data_list = input_data, y = outcome, cv_list = cv_list, p_metric = p_metric, algorithm = algorithm, n = iter)
    saveRDS(res,file.path("model_results", paste0(paste(outcome_name,integration,algorithm,p_metric,features,iter,sep = "_"),".rds")))
  } else {
    stop("Unkown integration! Choose either 'FFS' or 'ensemble'")
  }
} else {
  if (integration == "FFS"){
    res <- fit_forwardSelect(data_list = input_data, y = outcome, cv_list = cv_list, p_metric = p_metric, algorithm = algorithm, n = iter)
    saveRDS(res,file.path("model_results", paste0(paste(outcome_name,integration,algorithm,p_metric,features,iter,sep = "_"),".rds")))
  } else if (integration == "ensemble") {
    res <- fit_ensemble(data_list = input_data, y = outcome, cv_list = cv_list, p_metric = p_metric, algorithm = algorithm, n = iter)
    saveRDS(res,file.path("model_results", paste0(paste(outcome_name,integration,algorithm,p_metric,features,iter,sep = "_"),".rds")))
  } else {
    stop("Unkown integration! Choose either 'FFS' or 'ensemble'")
    }
}
