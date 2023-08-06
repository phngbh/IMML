#' Loading necessary data
#'
#' @description Loading all necessary data for building the package. The values
#' for the parameters can be TRUE or FALSE, if the respective data for the
#' modality should be loaded or not.
#'
#' @param clinical
#' @param metabolomics
#' @param methylomics
#' @param protemoics
#' @param transcriptomics
#'
#' @return Loads the selected data to the global environment.
#'
#' @author Ulrich Asemann

Load <- function(clinical = FALSE,
                 metabolomics = FALSE,
                 methylomics = FALSE,
                 protemoics = FALSE,
                 transcriptomics = FALSE,
                 genomics = FALSE) {
  # Function to load required data sheets and functions for building the package!
  # Only for use by builder of package

  # librarys
  library(devtools)
  library(roxygen2)
  library(readr)

  # remove exponential notation
  options(scipen = 999)#

  # General data used for everything
  sample_IDs <-
    readRDS("Add your path here!")
  phenotypes <-
    readRDS("Add your path here!")

  .GlobalEnv$sample_IDs <- sample_IDs
  .GlobalEnv$phenotypes <- phenotypes

  # DataPartitioning tables
  .GlobalEnv$training_IDs <-
    DataPartitioning(phenotypes, sample_IDs, seed = 124)
  .GlobalEnv$testing_IDs <-
    DataPartitioning(phenotypes, sample_IDs, type = "testing", seed = 124)
  .GlobalEnv$modeltraining_IDs <-
    DataPartitioning(phenotypes, sample_IDs, type = "modeltraining", seed = 124)


  # Clinical
  if (clinical) {
    clinical_processed <-
      readRDS("Add your path here!")

    # Transpose matrix, so its equal to the other data tables
    .GlobalEnv$clinical_processed <- clinical_processed
  }

  # Genomics
  if (genomics) {
    geno_annot <- read.csv(
      "Add your path here!",
      sep = "",
      check.names = F,
      header = F
    )
    geno_lowest_level_pathways <-
      readRDS("Add your path here!")

    .GlobalEnv$geno_annot <- geno_annot
    .GlobalEnv$geno_lowest_level_pathways <-
      geno_lowest_level_pathways
  }

  # Metabolomics
  if (metabolomics) {
    meta_lowest_level_pathways <-
      readRDS("Add your path here!")
    meta_anno_fil <-
      readRDS("Add your path here!")
    meta_pathway_list <-
      readRDS("Add your path here!")
    metabolomics_processed <-
      readRDS("Add your path here!")

    .GlobalEnv$meta_lowest_level_pathways <-
      meta_lowest_level_pathways
    .GlobalEnv$meta_pathway_list <- meta_pathway_list
    .GlobalEnv$meta_anno_fil <- meta_anno_fil
    .GlobalEnv$metabolomics_processed <- metabolomics_processed
  }

  # Methylomics
  if (methylomics) {
    meth_lowest_level_pathways <-
      readRDS("Add your path here!")
    methylomics_processed <-
      readRDS("Add your path here!")

    .GlobalEnv$meth_lowest_level_pathways <-
      meth_lowest_level_pathways
    .GlobalEnv$methylomics_processed <-
      methylomics_processed
  }

  # Protemoics
  if (protemoics) {
    prot_lowest_level_pathways <-
      readRDS("Add your path here!")
    prot_somamer_info_edited <-
      read.csv("Add your path here!")
    proteomics_processed <-
      readRDS("Add your path here!")

    .GlobalEnv$proteomics_processed <-
      proteomics_processed
    .GlobalEnv$prot_lowest_level_pathways <-
      prot_lowest_level_pathways
    .GlobalEnv$prot_somamer_info_edited <- prot_somamer_info_edited
  }

  # Transcriptomics
  if (transcriptomics) {
    trans_lowest_level_pathways <-
      readRDS("Add your path here!")
    trans_annotation_HumanHT_v3_final <-
      read.csv("Add your path here!")
    trans_gene_information <-
      readRDS("Add your path here!")
    transcriptomics_processed <-
      readRDS("Add your path here!")

    .GlobalEnv$trans_lowest_level_pathways <-
      trans_lowest_level_pathways
    .GlobalEnv$trans_annotation_HumanHT_v3_final <-
      trans_annotation_HumanHT_v3_final
    .GlobalEnv$trans_gene_information <- trans_gene_information
    .GlobalEnv$transcriptomics_processed <-
      transcriptomics_processed
  }


}
