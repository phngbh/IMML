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
  options(scipen = 999)


  # data to enviroment
  clinical_processed <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/R/IMLdata/data/Clinical_data/clinical_processed.rds"
    )
  sample_IDs <-
    readRDS("C:/Users/ulric/Desktop/Arbeit/R/IMLdata/supporting_files/sample_IDs.rds")
  phenotypes <-
    readRDS("C:/Users/ulric/Desktop/Arbeit/R/IMLdata/supporting_files/phenotypes.rds")
  transcriptomics_processed <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/R/IMLdata/data/Transcriptomics/transcriptomics_processed.rds"
    )
  metabolomics_processed <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/R/IMLdata/data/Metabolomics/metabolomics_processed.rds"
    )
  lowest_level_pathways <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/download/lowest_level_pathways.rds"
    )
  Annotation_HumanHT_v3_final <-
    read.csv(
      "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/download/Annotation_HumanHT-12v3_final.csv"
     )
  gene_information <-
    readRDS("C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/download/gene_IDs.rds")

  .GlobalEnv$clinical_processed <- clinical_processed
  .GlobalEnv$sample_IDs <- sample_IDs
  .GlobalEnv$phenotypes <- phenotypes
  .GlobalEnv$transcriptomics_processed <- transcriptomics_processed
  .GlobalEnv$metabolomics_processed <- metabolomics_processed
  .GlobalEnv$lowest_level_pathways <- lowest_level_pathways
  .GlobalEnv$annotation_HumanHT_v3_final <-
    Annotation_HumanHT_v3_final
  .GlobalEnv$gene_information <- gene_information

  # data_partitioning tables
  .GlobalEnv$training_IDs <-
    data_partitioning(phenotypes, sample_IDs)
  .GlobalEnv$testing_IDs <-
    data_partitioning(phenotypes, sample_IDs, type = "testing")

}
