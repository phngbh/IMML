#' load()
#'
#' @description Loading all necessary data, to build the package
#'
#' @return
#' @export
#'
#' @examples load()

load <- function(clinical = FALSE,
                 metabolomics = FALSE,
                 methylomics = FALSE,
                 protemoics = FALSE,
                 transcriptomics = FALSE) {
  # Function to load required data sheets and functions for building the package!
  # Only for use by builder of package

  # librarys
  library(devtools)
  library(roxygen2)
  library(readr)

  # remove exponential notation
  options(scipen = 999)


  # General data used for everything
  sample_IDs <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/supporting_files/sample_IDs.rds"
    )
  phenotypes <-
    readRDS(
      "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/supporting_files/phenotypes.rds"
    )

  .GlobalEnv$sample_IDs <- sample_IDs
  .GlobalEnv$phenotypes <- phenotypes

  # data_partitioning tables
  .GlobalEnv$training_IDs <-
    data_partitioning(phenotypes, sample_IDs)
  .GlobalEnv$testing_IDs <-
    data_partitioning(phenotypes, sample_IDs, type = "testing")


  # Clinical
  if (clinical) {
    clinical_processed <-
      readRDS(
        "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/data/Clinical_data/clinical_processed.rds"
      )

    # Transpose matrix, so its equa to the other data tables
    .GlobalEnv$clinical_processed <- t(clinical_processed)
  }

  # Genomics
  # coming soon!

  # Metabolomics
  if (metabolomics) {
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
    metabolomics_processed <-
      readRDS(
        "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/data/Metabolomics/metabolomics_processed.rds"
      )

    .GlobalEnv$meta_lowest_level_pathways <-
      meta_lowest_level_pathways
    .GlobalEnv$meta_pathway_list <- meta_pathway_list
    .GlobalEnv$meta_anno_fil <- meta_anno_fil
    .GlobalEnv$metabolomics_processed <- metabolomics_processed
  }

  # Methylomics
  if (methylomics) {
    meth_lowest_level_pathways <-
      readRDS(
        "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/download/fs_methylomics/meth_lowest_level_pathways.rds"
      )

    # methylomics_processed missing!

    .GlobalEnv$meth_lowest_level_pathways <-
      meth_lowest_level_pathways
  }

  # Protemoics
  if (protemoics) {
    prot_lowest_level_pathways <-
      readRDS(
        "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/download/fs_proteomics/prot_lowest_level_pathways.rds"
      )
    prot_somamer_info_edited <-
      read.csv(
        "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/download/fs_proteomics/prot_somamer_info_edited.csv"
      )

    # proteomics_processed missing!

    .GlobalEnv$prot_lowest_level_pathways <-
      prot_lowest_level_pathways
    .GlobalEnv$prot_somamer_info_edited <- prot_somamer_info_edited
  }

  # Transcriptomics
  if (transcriptomics) {
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
    transcriptomics_processed <-
      readRDS(
        "C:/Users/ulric/Desktop/Arbeit/Daten/IMLdata/data/Transcriptomics/transcriptomics_processed.rds"
      )

    .GlobalEnv$trans_lowest_level_pathways <-
      trans_lowest_level_pathways
    .GlobalEnv$trans_annotation_HumanHT_v3_final <-
      trans_annotation_HumanHT_v3_final
    .GlobalEnv$trans_gene_information <- trans_gene_information
    .GlobalEnv$transcriptomics_processed <-
      transcriptomics_processed
  }


}
