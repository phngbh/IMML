#' Feature Selection Clinical Data
#'
#' @description Feature selection for clinical data using elastic net models.
#'
#' @param trainIDs Set of training IDs from the `DataPartitioning()` function
#'   for the feature selection.
#' @param testIDs Set of testing IDs from the `DataPartitioning()` function for
#'   the feature selection.
#' @param dataIDs A data.frame with samples as rows and the data modalities as
#'   columns. Row names represent sample IDs. It holds the data IDs of a sample
#'   for each modality. If for a sample there is no data for a modality, it has
#'   to be indicated by NA.
#' @param phenotypeIDs A data.frame with samples as rows and sample IDs as row
#'   names. Columns are phenotypes of interest. If for a sample no information
#'   about a phenotype is available, it has to be indicated by NA.
#' @param phenotype The name of the column in phenotypeIDs, which will be used
#'   in the analysis.
#' @param clinicalData A data.frame with samples as rows and sample IDs as row
#'   names. Columns are the clinical variables and are named accordingly.
#' @param seed The seed used for random number generation. Using the same seed
#'   ensures reproducibility.
#'
#' @return Returns clincalData subset to the selected features.
#'
#' @author Ulrich Asemann & Wilhelm Glaas
#'
#' @export

FsClinical <-
  function(trainIDs,
           testIDs,
           dataIDs,
           phenotypeIDs,
           phenotype,
           clinicalData,
           seed = 123) {

    # Save row names as a column
    dataIDs <- dataIDs %>%
      tibble::rownames_to_column(var = "sampleIDs")
    phenotypeIDs <- phenotypeIDs %>%
      tibble::rownames_to_column(var = "sampleIDs")

    merged <- merge(dataIDs, phenotypeIDs, by = "sampleIDs") %>%
      select(Clinical, any_of(phenotype)) %>% drop_na()

    indicator <- as.vector(merged)[[2]]
    names(indicator) <- as.vector(merged)[[1]]
    indicator <-
      ifelse(indicator == 1, "One", "Zero") %>%
      factor(levels = c("One", "Zero"))

    # Do elastic net
    # Build control
    myControl <- trainControl(
      method = "repeatedcv",
      number = 5,
      repeats = 4,
      savePredictions = "final",
      classProbs = TRUE,
      summaryFunction = CaretLogLoss,
      sampling = NULL,
      allowParallel = TRUE
    )

    # Create variables
    perfList = list()
    var = list()

    # Loop through all partitions
    for (i in 1:length(trainIDs)) {
      cat("Fold", i, "\n")

      # Select the IDs for the partition
      trainID <- unlist(trainIDs[[i]])
      testID <- unlist(testIDs[[i]])

      # Using the rownames of the selected IDs for x also for y
      # NA values are no longer used this way
      xTrain <-
        clinicalData[as.character(trainID),] %>% drop_na()

      yTrain <- indicator[as.character(rownames(xTrain))]

      xTest <- clinicalData[as.character(testID),] %>% drop_na()

      yTest <- indicator[rownames(xTest)]

      weights <-
        ifelse(yTrain == "One", table(yTrain)[[2]] / table(yTrain)[[1]], 1)

      # Building the fit
      set.seed(seed)
      fit <- caret::train(
        x = xTrain,
        y = yTrain,
        method = "glmnet",
        metric = "MyLogLoss",
        tuneLength = 20,
        weights = weights,
        maximize = FALSE,
        trControl = myControl,
        importance = TRUE
      )

      # Building the prediction
      pred <-
        predict(fit, xTest, s = "lambda.min", type = "prob")$One

      roc <-
        roc(
          response = yTest,
          predictor = pred,
          levels = c("Zero", "One")
        )
      auc <- auc(roc)
      perfList[[i]] <- auc
      varImp <- varImp(fit)$importance
      vari <- varImp$Overall
      names(vari) <- rownames(varImp)
      vari <- vari[order(vari, decreasing = T)]

      var[[i]] <- vari
    }


    # Select best features
    varFil <- lapply(var, function(x)
      names(x[x > 0]))
    varRra <- aggregateRanks(varFil)
    varRra$adjP <- varRra$Score * 100
    varRra$adjP <- p.adjust(varRra$adjP, "fdr")
    varSel <- varRra %>% dplyr::filter(adjP < 0.05)
    varSel <- unname(unlist(select(varSel, "Name")))

    # Save the data
    #saveRDS(dplyr::select(clinicalData, all_of(varSel)),
    #        "clinical_selected.rds")

    result <- select(clinicalData, all_of(varSel))

    return(result)
  }
