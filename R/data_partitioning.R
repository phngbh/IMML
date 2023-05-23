#' Data Partitioning
#'
#' @description Returns a list with all training IDs and
#' a set of shuffled IDs for feature selection.
#'
#' @param phenotype_IDs A table holding the clinical IDs with a variable indicating,
#' if the disease occurred or not.
#' @param data_IDs A table holding all the clinical IDs with the respective IDs
#' for each modality.
#' @param partitioning A value, that signals what percentage of IDs are used for
#' training in the feature selection.
#' @param type This variable signals, if the training or testing set of the
#' feature selection is returned.
#' @param amount The number of partitions that are build for feature selection
#' @param seed The possibility to change the seed, for different results in the
#' part of the ID shuffling for feature selection.
#'
#' @return Returns a list, holding all IDs for model training and a set of
#' shuffled IDs depending on the chosen amount for each modality. The different
#' sets for the modalities are hold separately in extra sublists.
#'
#' @author Ulrich Asemann

data_partitioning <-
  function(phenotype_IDs,
           data_IDs,
           partitioning = 0.8,
           amount = 100,
           type = "training",
           seed = 123) {
    # Set seed, achieving the same results and the complementary testing IDs
    # for each partition
    set.seed(seed)

    # Defining variables
    tmp_list <- list()
    training_IDs <- list()
    single_feature_list <- list()
    feature_list <- list()
    return_list <- list()

    # List of feature selection modalities
    list_names <- colnames(data_IDs)

    # Selecting IDs used for training and feature selection
    IDs <-
      phenotype_IDs %>% dplyr::select(inc3) %>% drop_na() %>% rownames() %>% as.double()

    # Data frame with the selected IDs
    new_sample_IDs <- data_IDs %>% filter(Clinical %in% IDs)


    # Building the modeltraining IDs list
    # Selecting clinical IDs with IDs for each modality
    training_sample_IDs <-
      new_sample_IDs %>% drop_na() %>% dplyr::select(Clinical)

    # Building a list for the modeltraining IDs
    for (i in 1:nrow(training_sample_IDs)) {
      training_IDs <- c(training_IDs, training_sample_IDs[i, ])
      tmp_name <- paste("Modeltraining", i)
      names(training_IDs)[i] <- tmp_name
    }

    # Remove the modeltraining IDs from the new_sample_IDs
    remove_IDs <- unlist(training_IDs)
    new_sample_IDs <-
      new_sample_IDs %>% filter(!(Clinical %in% remove_IDs))

    # Building a feature selection list for each modality
    for (i in 1:length(list_names)) {
      # Select one feature and remove all NA values from the new_sample_IDs list
      frame <-
        new_sample_IDs %>% dplyr::select(list_names[i], Clinical) %>% drop_na() %>% dplyr::select(Clinical)

      # Creating data partitions with "caret"
      data_part <-
        caret::createDataPartition(unlist(frame), times = amount, p = partitioning)

      # Go through all 100 partitions
      for (j in 1:length(data_part)) {
        # Build the feature selection lists, depending on training or testing
        if (type == "testing") {
          tmp_frame <- frame[-(unlist(data_part[j])), ]
        } else if (type == "training") {
          tmp_frame <- frame[(unlist(data_part[j])), ]
        } else {
          return("Please set the type to testing or training!")
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

    # Put the modeltraining IDs list and feature selection lists into one
    if (type == "testing") {
      return_list <-
        list("Modeltraining IDs" = training_IDs,
             "Testing Feature Selection IDs" = feature_list)
    } else if (type == "training") {
      return_list <-
        list("Modeltraining IDs" = training_IDs,
             "Training Feature Selection IDs" = feature_list)
    }

    # Return the completed list
    return(return_list)
  }
