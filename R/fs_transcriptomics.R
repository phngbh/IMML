#' Feature Selection for Transcriptomics
#'
#' @param IDs
#' @param data_IDs
#' @param transcriptomics_data
#'
#' @return
#' @export
#'
#' @examples
fs_transcriptomics <- function(IDs, data_IDs, transcriptomics_data) {
  #
  # Getting the ID sets from the transcriptomics
  transcriptomics_IDs <- IDs$`Feature Selection IDs`$Transcriptomics

  # Going through each set of IDs
  for (i in 1:length(transcriptomics_IDs)) {
    # Selecting the data with one set of the IDs
    tmp_IDs <- unlist(transcriptomics_IDs[[i]])

    tmp_transcriptomics_IDs <-
      data_IDs %>% filter(Clinical %in% tmp_IDs) %>%
      select(Transcriptomics) %>% as.list() %>% unlist()

    tmp_transcriptomics <-
      transcriptomics_data %>% as.data.frame() %>%
      select(suppressWarnings(one_of(
        as.character(tmp_transcriptomics_IDs)
      )))

    # Sort the column names (will be removed eventually)
    tmp_transcriptomics <-
      tmp_transcriptomics[, order(names(tmp_transcriptomics))]


    return(tmp_transcriptomics)
  }

  return(length(transcriptomics_IDs))

}
