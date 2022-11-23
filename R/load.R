#' load()
#'
#' @description Loading all necessary data, to build the package
#'
#' @return
#' @export
#'
#' @examples load()

load <- function() {
  # Function to load required data sheets and functions for building the package!
  # Only for use by builder of package

  # librarys
  library(devtools)
  library(roxygen2)
  library(readr)

  # remove exponential notation
  # options(scipen = 999)


  # General data used for everything
  clinical_processed <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/data/Clinical_data/clinical_processed.rds"
    )
  sample_IDs <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/supporting_files/sample_IDs.rds"
    )
  phenotypes <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/supporting_files/phenotypes.rds"
    )
  transcriptomics_processed <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/data/Transcriptomics/transcriptomics_processed.rds"
    )
  metabolomics_processed <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/data/Metabolomics/metabolomics_processed.rds"
    )

  # Feature Selection Transcriptomics
  trans_lowest_level_pathways <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/download/fs_transcriptomics/trans_lowest_level_pathways.rds"
    )
  trans_annotation_HumanHT_v3_final <-
    read.csv(
      "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/download/fs_transcriptomics/trans_annotation_HumanHT-12v3_final.csv"
    )
  trans_gene_information <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/download/fs_transcriptomics/trans_gene_IDs.rds"
    )

  # Feature Selection Metabolomics
  meta_lowest_level_pathways <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/download/fs_metabolomics/meta_lowest_level_pathways.rds"
    )
  meta_anno_fil <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/download/fs_metabolomics/meta_anno_fil.rds"
    )
  meta_pathway_list <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/download/fs_metabolomics/meta_pathway_list.rds"
    )



  # General data
  .GlobalEnv$clinical_processed <- clinical_processed
  .GlobalEnv$sample_IDs <- sample_IDs
  .GlobalEnv$phenotypes <- phenotypes
  .GlobalEnv$transcriptomics_processed <- transcriptomics_processed
  .GlobalEnv$metabolomics_processed <- metabolomics_processed

  # Transcriptomics data
  .GlobalEnv$trans_lowest_level_pathways <-
    trans_lowest_level_pathways
  .GlobalEnv$trans_annotation_HumanHT_v3_final <-
    trans_annotation_HumanHT_v3_final
  .GlobalEnv$trans_gene_information <- trans_gene_information

  # Metabolomics data
  .GlobalEnv$meta_lowest_level_pathways <- meta_lowest_level_pathways
  .GlobalEnv$meta_pathway_list <- meta_pathway_list
  .GlobalEnv$meta_anno_fil <- meta_anno_fil

  # data_partitioning tables
  .GlobalEnv$training_IDs <-
    data_partitioning(phenotypes, sample_IDs)
  .GlobalEnv$testing_IDs <-
    data_partitioning(phenotypes, sample_IDs, type = "testing")

}
