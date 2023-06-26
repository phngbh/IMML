#' Data Partitioning
#'
#' @description Returns a list with all training IDs and
#' a set of shuffled IDs for feature selection.
#'
#' @param phenotypeIDs A table holding the clinical IDs with a variable indicating,
#' if the disease occurred or not.
#' @param dataIDs A table holding all the clinical IDs with the respective IDs
#' for each modality.
#' @param partitioning A value, that signals what percentage of IDs are used for
#' training in the feature selection.
#' @param type This variable signals, if the training or testing set of the
#' feature selection is returned or the modeltraining sets.
#' @param numPartitions The number of partitions that are build for feature selection
#' @param seed The possibility to change the seed, for different results in the
#' part of the ID shuffling for feature selection.
#'
#' @return Returns a list, holding all IDs for model training and a set of
#' shuffled IDs depending on the chosen amount of partitions for each modality.
#' The different sets for the modalities are hold separately in extra sublists.
#'
#' @author Ulrich Asemann

DataPartitioning <-
  function(phenotypeIDs,
           dataIDs,
           partitioning = 0.8,
           numPartitions = 100,
           type = "training",
           seed = 123) {
    # Set seed, to achieve same results
    set.seed(seed)

    # Defining variables
    tmpList <- list()
    tmpList2 <- list()
    trainingIDs <- list()
    singleFeatureList <- list()
    featureList <- list()
    returnList <- list()

    # List of feature selection modalities
    listNames <- colnames(dataIDs)

    # Selecting IDs used for training and feature selection
    IDs <-
      phenotypeIDs %>% dplyr::select(inc3) %>% drop_na() %>% rownames() %>% as.double()

    # Data frame with the selected IDs
    newSampleIDs <- dataIDs %>% filter(Clinical %in% IDs)


    # Building the model training IDs list
    # Selecting clinical IDs with IDs for each modality
    trainingSampleIDs <-
      newSampleIDs %>% drop_na() %>% dplyr::select(Clinical)


    # Building a list with the modeltraining IDs
    for (i in 1:nrow(trainingSampleIDs)) {
      trainingIDs <- c(trainingIDs, trainingSampleIDs[i,])
      # tmpName <- paste("Modeltraining", i)
      # names(trainingIDs)[i] <- tmpName
    }

    if (type == "modeltraining") {
      # Create partitions for the modeltraining
      dataPartitions <-
        caret::createDataPartition(unlist(trainingSampleIDs),
                                   times = 100,
                                   p = 0.8)

      # Create the partitions with the clinical IDs
      for (i in 1:length(dataPartitions)) {
        # Select the IDs
        tmpList <- trainingIDs[unlist(dataPartitions[i])]

        # Build the lists holding the IDs
        partitionTrainingIDs <-
          split(tmpList, cut(seq_along(tmpList), 4, labels = FALSE))
        partitionTestingIDs <-
          trainingIDs[-unlist(dataPartitions[i])]

        # Putting the lists together
        tmpList <-
          list(c(
            "Fold" = partitionTrainingIDs,
            list("Fold.5" = partitionTestingIDs)
          ))

        returnList <- c(returnList, tmpList)
        names(returnList)[i] <- paste("Modeltraining", i)

      }
      # Return the result
      returnList <- list("Modeltraining" = returnList)

      return(returnList)
    }

    # Remove the modeltraining IDs from the newSampleIDs
    removeIDs <- unlist(trainingIDs)
    newSampleIDs <-
      newSampleIDs %>% filter(!(Clinical %in% removeIDs))

    # Building a feature selection list for each modality
    for (i in 1:length(listNames)) {
      # Select one feature and remove all NA values from the newSampleIDs list
      frame <-
        newSampleIDs %>% dplyr::select(listNames[i], Clinical) %>% drop_na() %>% dplyr::select(Clinical)

      # Creating data partitions with "caret"
      dataPartitions <-
        caret::createDataPartition(unlist(frame), times = numPartitions, p = partitioning)

      # Go through all 100 partitions
      for (j in 1:length(dataPartitions)) {
        # Build the feature selection lists, depending on training or testing
        if (type == "testing") {
          tmpFrame <- frame[-(unlist(dataPartitions[j])),]
        } else if (type == "training") {
          tmpFrame <- frame[(unlist(dataPartitions[j])),]
        } else {
          return("Please set the type to testing or training!")
        }

        # Building the list of the current partition
        for (k in 1:nrow(tmpFrame)) {
          tmpList <- c(tmpList, tmpFrame[k,])
          newName <- paste(listNames[i], j, k)
          names(tmpList)[k] <- newName
        }

        # Adding the tmpList to the singleFeatureList
        singleFeatureList[[j]] <- tmpList
        names(singleFeatureList)[j] <- paste(listNames[i], j)

        # Empty the tmpList
        tmpList <- list()

      }

      # Add the singleFeatureList to their respective feature
      assign(listNames[i], singleFeatureList)
    }

    # Put all singleFeatureLists into one list
    for (i in 1:length(listNames)) {
      featureList[[i]] <- get(listNames[i])
      names(featureList)[i] <- listNames[i]
    }

    # Put the modeltraining IDs list and feature selection lists into one
    if (type == "testing") {
      returnList <-
        list("Testing Feature Selection IDs" = featureList)
    } else if (type == "training") {
      returnList <-
        list("Training Feature Selection IDs" = featureList)
    } else if (type == "modeltraining") {
      returnList <- list("Modeltraining IDs" = trainingIDs)
    }

    # Return the completed list
    return(returnList)
  }
