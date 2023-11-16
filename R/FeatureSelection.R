#' Feature Selection Superfunction
#'
#' @description Performs feature selection on the modalities specified by the
#'   user.
#'
#' @param dataPart Partitioned data output by the `DataPartitioning()` function.
#' @param modalities A vector specifying the modalities, which are to be used in
#'   feature selection. E.g.: `c('Clinical', 'Genomics', 'Proteomics')`.
#'   Accepted modalities are Genomics, Transcriptomics, Proteomics,
#'   Metabolomics, Methylomics and Clinical.
#' @param dataList A named list of the data.frames holding the respective data
#'   for each modality. List elements are named after the repsective modality.
#'   See @param modalities for accepted names.
#' @param dataIDs A data.frame with samples as rows and the data modalities as
#'   columns. Row names represent sample IDs. It holds the data IDs of a sample
#'   for each modality. If for a sample there is no data for a modality, it has
#'   to be indicated by NA.
#' @param phenotypeIDs A data.frame with samples as rows and sample IDs as row
#'   names. Columns are phenotypes of interest. If for a sample no information
#'   about a phenotype is available, it has to be indicated by NA.
#' @param phenotype The name of the column in phenotypeIDs, which will be used
#'   in the analysis.
#' @param pathways Gene pathways used in GSEA. !!!add more details!!!
#' @param seed The seed used for random number generation. Using the same seed
#'   ensures reproducibility.
#'
#' @return The result of another function, depending on the chosen modality
#'
#' @author Ulrich Asemann & Wilhelm Glaas
#'
#' @export

FeatureSelection <-
  function(dataPart,
           modalities,
           dataList,
           dataIDs,
           phenotypeIDs,
           phenotype,
           seed = 123) {
    #Checking for correct input of modalities
    acceptedMods = c('Genomics', 'Transcriptomics', 'Proteomics',
                     'Metabolomics', 'Methylomics', 'Olink', 'Clinical')
    for (mod in modalities) {
      if(!(mod %in% names(dataPart$`Testing Feature Selection IDs`))) {
        stop(paste0(mod, ' is not present in the data partition!'))
      }
      if(!(mod %in% acceptedMods)) {
        stop(paste0(mod, ' is not an accepted modality!'))
      }
    }

    fsTrain <- dataPart$trainingFeatureSelectionIDs
    fsTest <- dataPart$testingFeatureSelectionIDs


    returnList <- case_when(modalities == 'Genomics' ~
                              FsGenomics(),
                            modalities == 'Transcriptomics' ~
                              FsTranscriptomics(trainIDs = fsTrain$Transciptomics,
                                                testIDs = fsTest$Transcriptomics,
                                                ),
                            modalities == 'Proteomics' ~
                              FsProteomics(),
                            modalities == 'Metabolomics' ~
                              FsMetabolomics(),
                            modalities == 'Methylomics' ~
                              FsMethylomics(),
                            modalities == 'Clinical' ~
                              FsClinical(trainIDs = fsTrain$Clinical,
                                         testIDs = fsTest$Clinical,
                                         dataIDs = dataIDs,
                                         phenotypeIDs = phenotypeIDs,
                                         phenotype = phenotype,
                                         clinicalData = dataList$`Clinical`,
                                         seed = seed)
    )


    returnList = as.list(returnList)
    names(returnList) = modalities

    return(returnList)
}
