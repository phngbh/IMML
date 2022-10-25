#' Title
#'
#' @param clinical_IDs
#' @param data
#'
#' @return
#' @export
#'
#' @examples


data_partitioning <- function(clinical_IDs, data) {
  # Variables
  training_IDs <- list()
  feature_list <- list()
  return_list <- list()
  # List of feature selection parts
  list_names <-
    c("Genomics",
      "Transcriptomics",
      "Proteomics",
      "Metabolomics",
      "Methylomics")

  # The IDs used for training and feature selection
  IDs <-
    clinical_IDs %>% select(inc3) %>% drop_na() %>% rownames() %>% as.double()

  # A data frame only with clinical_IDs, that are used for training and feature selection
  new_sample_IDs <- data %>% filter(Clinical %in% IDs)


  # Building a list with the training IDs

  # Only using IDs, that have values in all columns and select the training IDs
  training_sample_IDs <-
    new_sample_IDs %>% drop_na() %>% select(Clinical)

  # Build a list with the 100 training IDs
  for (i in 1:100) {
    training_IDs <- c(training_IDs, training_sample_IDs[i,])
    tmp_name <- paste("Training", i)
    names(training_IDs)[i] <- tmp_name
  }

  # Remove the training_IDs from the new_sample_IDs
  remove_IDs <- unlist(training_IDs)
  new_sample_IDs <-
    new_sample_IDs %>% filter(!(Clinical %in% remove_IDs))

  # Building a feature selection list for each single feature
  for (i in 1:length(list_names)) {
    # Select one feature and remove all NA values from the new_sample_IDs list
    frame <-
      new_sample_IDs %>% select(list_names[i], Clinical) %>% drop_na() %>% select(Clinical)

    tmp_list <- list()
    tmp_name <- list_names[i]

    for (j in 1:nrow(frame)) {
      tmp_list <- c(tmp_list, frame[j,])
      new_name <- paste(tmp_name, j)
      names(tmp_list)[j] <- new_name

    }

    assign(list_names[i], tmp_list)
  }

  # Put all feature selection ID lists into one list
  for (i in 1:length(list_names)) {
    feature_list[[i]] <- get(list_names[i])
    names(feature_list)[i] <- list_names[i]
  }

  # Put training IDs list and feature selection IDs lists into one
  return_list <-
    list("Training IDs" = training_IDs,
         "Feature Selection IDs" = feature_list)

  return(return_list)
}
