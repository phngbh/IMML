#' Data Partitioning
#'
#' @description Builds a list data partitions for the provided sets of IDs.
#'
#' @param phenotypeIDs A data.frame with samples as rows and sample IDs as row
#'   names. Columns are phenotypes of interest. If for a sample no information
#'   about a phenotype is available, it has to be indicated by NA.
#' @param dataIDs A data.frame with samples as rows and the data modalities as
#'   columns. Row names represent sample IDs. It holds the data IDs of a sample
#'   for each modality. If for a sample there is no data for a modality, it has
#'   to be indicated by NA.
#' @param phenotype The name of the column in phenotypeIDs, which will be used
#'   in the analysis.
#' @param partitioning A value, that determines what percentage of samples are
#'   used for training during feature selection and model training.
#' @param numPartitions The number of partitions that are build for feature
#'   selection and model training.
#' @param k The amount of folds used for k-fold cross validation during model
#'   training.
#' @param iter The amount of iterations for k-fold cross validation during model
#'   training.
#' @param seed The seed used for random number generation. Using the same seed
#'   ensures reproducibility.
#'
#' @return The function return three lists, each containing partitions of the
#'   sample IDs for feature selection testing/training and model training
#'   respectively.
#'
#' @author Ulrich Asemann & Wilhelm Glaas

DataPartitioning <-
  function(phenotypeIDs,
           dataIDs,
           phenotype,
           partitioning = 0.8,
           numPartitions = 100,
           k = 5,
           iter = 100,
           seed = 123) {
    # Setting seed
    set.seed(seed)

    # Defining variables
    trainingIDs <- list()
    trainSingleFeatureList <- list()
    testSingleFeatureList <- list()
    testFeatureList <- list()
    trainFeatureList <- list()
    finalList <- list()

    # List of feature selection modalities
    listNames <- colnames(dataIDs)

    # Selecting IDs used for training and feature selection
    IDs <-
      phenotypeIDs %>%
        dplyr::select(any_of(phenotype)) %>%
        drop_na() %>%
        rownames() %>%
        as.double()

    # Save row names as a column
    dataIDs <- dataIDs %>% tibble::rownames_to_column(var = "sampleID")

    # Data frame with the selected IDs
    newSampleIDs <- dataIDs %>% filter(sampleID %in% IDs)


    # Building the model training IDs list

    # Selecting sample IDs with IDs for each modality
    trainingSampleIDs <-
      newSampleIDs %>% drop_na() %>% dplyr::select(sampleID)


    # Building a list with the modeltraining IDs
    trainingIDs <- as.double(trainingSampleIDs$sampleID)

    iterations <-
      caret::createDataPartition(trainingIDs, times = iter, p = partitioning)

    returnList <- list()
    for (i in 1:iter) {
      tmpFoldList <- list()

      testIDs <- trainingIDs[-iterations[[i]]]
      trainIDs <- createFolds(trainingIDs[iterations[[i]]], k = k)
      trainIDs <- lapply(trainIDs, function(x) trainingIDs[x])

      tmpFoldList[['Training']] <- trainIDs
      tmpFoldList[['Testing']] <- testIDs

      returnList[[paste0('Iteration', i)]] <- tmpFoldList
    }

    # append to the final list
    finalList <- c(finalList, "Modeltraining" = list(returnList))


    # Remove the modeltraining IDs from the newSampleIDs
    newSampleIDs <-
      newSampleIDs %>% filter(!(sampleID %in% trainingIDs))


    # Building a feature selection list for each modality
    for (i in 1:length(listNames)) {
      # Select one feature and remove all NA values from the newSampleIDs list
      frame <-
        newSampleIDs %>% dplyr::select(listNames[i]) %>% drop_na()

      # Creating data partitions with "caret"
      dataPartitions <-
        caret::createDataPartition(unlist(frame),
                                   times = numPartitions, p = partitioning)

      # Iterating over the partitions
      for (j in 1:length(dataPartitions)) {

        # Build the feature selection lists
        testTmpList <- unlist(frame[-dataPartitions[[j]], ])
        trainTmpList <- unlist(frame[dataPartitions[[j]], ])

        # Adding the tmpLists to the singleFeatureLists
        testSingleFeatureList[[j]] <- unname(testTmpList)
        names(testSingleFeatureList)[j] <- paste0(listNames[i], j)

        trainSingleFeatureList[[j]] <- unname(trainTmpList)
        names(trainSingleFeatureList)[j] <- paste0(listNames[i], j)

      }

      # Add the testSingleFeatureLists to their respective featureList
      testFeatureList[i] <- list(testSingleFeatureList)
      names(testFeatureList)[i] <- listNames[i]

      trainFeatureList[i] <- list(trainSingleFeatureList)
      names(trainFeatureList)[i] <- listNames[i]

    }

    finalList <- c(finalList,
                   "Testing Feature Selection IDs" = list(testFeatureList),
                   "Training Feature Selection IDs" = list(trainFeatureList))

    # Return the completed list
    return(finalList)
  }
