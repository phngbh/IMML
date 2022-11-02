#' Data Partitioning
#'
#' @description Returns a list with all training IDs and
#' a set of 100 shuffled feature selection IDs.
#'
#' @author Ulrich Asemann
#'
#' @param phenotype_IDs
#' @param data_IDs
#' @param partitioning
#' @param type
#'
#' @return
#' @export
#'
#' @examples
data_partitioning <-
  function(phenotype_IDs,
           data_IDs,
           partitioning = 0.8,
           type = "training",
           seed = 123) {
    #
    # Set seed, to get the same results for training and testing
    # Can be modified by the user
    set.seed(seed)

    # Variables
    tmp_list <- list()
    training_IDs <- list()
    single_feature_list <- list()
    feature_list <- list()
    return_list <- list()

    # List of feature selection parts
    list_names <- colnames(data_IDs)

    # The IDs used for training and feature selection
    IDs <-
      phenotype_IDs %>% select(inc3) %>% drop_na() %>% rownames() %>% as.double()

    # A data frame only with phenotype_IDs, that are used for training and feature selection
    new_sample_IDs <- data_IDs %>% filter(Clinical %in% IDs)


    # Building a list with the training IDs

    # Only using IDs, that have values in all columns and select the training IDs
    training_sample_IDs <-
      new_sample_IDs %>% drop_na() %>% select(Clinical)

    # Build a list with the 100 training IDs
    for (i in 1:nrow(training_sample_IDs)) {
      training_IDs <- c(training_IDs, training_sample_IDs[i, ])
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

      # Create a data partitioning with "caret"
      data_part <-
        createDataPartition(unlist(frame), times = 100, p = partitioning)

      # Go through each of the 100 partitions
      for (j in 1:length(data_part)) {
        tmp_frame <- frame[(unlist(data_part[j])), ]


        # In case of testing, remove the training data from the partition
        if (type == "testing") {
          tmp_frame <- frame %>% slice(-(unlist(data_part[j])))
        }

        # Building the list of the current partition
        for (k in 1:nrow(tmp_frame)) {
          tmp_list <- c(tmp_list, tmp_frame[k, ])
          new_name <- paste(list_names[i], j, k)
          names(tmp_list)[k] <- new_name
        }

        # Adding the tmp_list to the single_feature_list
        single_feature_list[[j]] <- tmp_list
        names(single_feature_list)[j] <- paste(list_names[i], j)

        # Empty the tmp_list
        tmp_list <- list()

      }

      # Add the single_feature_list to their feature
      assign(list_names[i], single_feature_list)
    }

    # Put all single_feature_lists into one list
    for (i in 1:length(list_names)) {
      feature_list[[i]] <- get(list_names[i])
      names(feature_list)[i] <- list_names[i]
    }

    # Put training_IDs and feature_list into one
    return_list <-
      list("Training IDs" = training_IDs,
           "Feature Selection IDs" = feature_list)

    # Return the completed list
    return(return_list)
  }
