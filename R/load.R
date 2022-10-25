#' load()
#'
#' @description Loading all necessary data, to build the package
#'
#' @return
#' @export
#'
#' @examples load()

load <- function(){

  # Function to load required data sheets and functions for building the package!
  # Only for use by builder of package

  clinical_processed <- readRDS("C:/Users/ulric/Desktop/Arbeit/R/IMLdata/data/Clinical_data/clinical_processed.rds")
  sample_IDs <- readRDS("C:/Users/ulric/Desktop/Arbeit/R/IMLdata/supporting_files/sample_IDs.rds")
  phenotypes <- readRDS("C:/Users/ulric/Desktop/Arbeit/R/IMLdata/supporting_files/phenotypes.rds")
  transcriptomics_processed <- readRDS("C:/Users/ulric/Desktop/Arbeit/R/IMLdata/data/Transcriptomics/transcriptomics_processed.rds")
  metabolomics_processed <- readRDS("C:/Users/ulric/Desktop/Arbeit/R/IMLdata/data/Metabolomics/metabolomics_processed.rds")

  .GlobalEnv$clinical_processed <- clinical_processed
  .GlobalEnv$sample_IDs <- sample_IDs
  .GlobalEnv$phenotypes <- phenotypes
  .GlobalEnv$transcriptomics_processed <- transcriptomics_processed
  .GlobalEnv$metabolomics_processed <- metabolomics_processed

  library(devtools)
  library(roxygen2)
}
