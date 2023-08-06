#' Feature Selection Clinical Data
#'
#' @description The function for the feature selection for the clinical data.
#'
#' @param trainIDs Set of training IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param testIDs Set of testing IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param dataIDs A table holding all the clinical IDs with the respective IDs
#' for each modality.
#' @param phenotypeIDs A table holding the clinical IDs with a variable indicating,
#' if the disease occurred or not.
#' @param clinicalData A table holding the clinical data.
#'
#' @return Returns the clinical data table only with the significant selected features.
#'
#' @author Ulrich Asemann

FsClinical <-
  function(trainIDs,
           testIDs,
           dataIDs,
           phenotypeIDs,
           clinicalData) {
    # Getting IDs
    trainClinicalIDs <-
      trainIDs$`Training Feature Selection IDs`$Clinical
    testClinicalIDs <-
      testIDs$`Testing Feature Selection IDs`$Clinical

    filteredPhenotypeIDs <- filter(phenotypeIDs, !is.na(inc3))

    indicator <- filteredPhenotypeIDs$inc3

    names(indicator) <- rownames(filteredPhenotypeIDs)

    indicator <- na.omit(indicator)

    indicator <-
      ifelse(indicator == 1, "One", "Zero") %>% factor(levels = c("One", "Zero"))

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
    for (i in 1:length(trainClinicalIDs)) {
      cat("Fold", i, "\n")

      # Select the IDs for the partition
      trainID <- unlist(trainClinicalIDs[[i]])
      testID <- unlist(testClinicalIDs[[i]])

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
      set.seed(993)
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

    # Set seed
    set.seed(993)

    # Select best features
    varFil <- lapply(var, function(x)
      names(x[x > 0]))
    varRra <- aggregateRanks(varFil)
    varRra$adjP <- varRra$Score * 100
    varRra$adjP <- p.adjust(varRra$adjP, "fdr")
    varSel <- varRra %>% dplyr::filter(adjP < 0.05)
    varSel <- as.numeric(unname(unlist(select(varSel, "Name"))))

    # Save the data
    saveRDS(dplyr::select(clinicalData, all_of(varSel)),
            "clinical_selected.rds")

    return("Elasticnet clinical done!")
  }
