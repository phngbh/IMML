#' SelectModeltrainingIDs
#'
#' @description
#'
#' @param dataFrame
#' @param modeltrainingIDs
#' @param dataIDs
#' @param modality
#'
#' @return the modality data frame with only the modeltraining IDs
#'
#' @author Ulrich Asemann
SelectModeltrainingIDs <-
  function(dataFrame,
           modeltrainingIDs,
           dataIDs,
           modality) {
    # Check if the chosen modality is right
    if (!(
      modality %in% c(
        'Clinical',
        'Genomics',
        'Metabolomics',
        'Methylomics',
        'Proteomics',
        'Transcriptomics'
      )
    )) {
      return(print("The selected modality is wrong!"))
    }

    # Turn the data into a data frame
    dataFrame <- as.data.frame(dataFrame)

    # Get the modeltraining IDs
    unlistedIDs <- unlist(modeltrainingIDs) %>% unique()

    # Turn the cinical IDs into the modality IDs
    modalityIDs <-
      dataIDs[match(unlistedIDs, dataIDs$Clinical), ] %>% dplyr::select(modality) %>% unlist()

    # Create the data frame with only the modeltraining IDs
    selectedDataFrame <-
      dataFrame[, !is.na(match(colnames(dataFrame), modalityIDs))]


    return(selectedDataFrame)
  }
