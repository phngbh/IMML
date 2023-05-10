#' Feature Selection Clinical Data
#'
#' @description The function for the feature selection for the clinical data.
#'
#' @param train_IDs Set of training IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param test_IDs Set of testing IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param data_IDs A table holding all the clinical IDs with the respective IDs
#' for each modality.
#' @param phenotype_IDs A table holding the clinical IDs with a variable indicating,
#' if the disease occurred or not.
#' @param clinical_data A table holding the clinical data.
#'
#' @return Returns the clinical data table only with the significant selected features.
#'
#' @author Ulrich Asemann

fs_clinical <-
  function(train_IDs,
           test_IDs,
           data_IDs,
           phenotype_IDs,
           clinical_data) {
    #
    # Getting IDs
    train_clinical_IDs <-
      train_IDs$`Feature Selection IDs`$Clinical
    test_clinical_IDs <-
      test_IDs$`Feature Selection IDs`$Clinical

    mnsi <- filter(phenotype_IDs,!is.na(inc3))

    inc <- mnsi$inc3

    names(inc) <- rownames(mnsi)

    inc <- na.omit(inc)

    inc <-
      ifelse(inc == 1, "One", "Zero") %>% factor(levels = c("One", "Zero"))

    # Do elastic net
    # Build control
    my_control <- trainControl(
      method = "repeatedcv",
      number = 5,
      repeats = 4,
      savePredictions = "final",
      classProbs = TRUE,
      summaryFunction = caretLogLoss,
      sampling = NULL,
      allowParallel = T
    )

    # return(my_control)

    perf_list = list()
    var = list()


    # length(resamples)
    for (i in 1:100) {
      cat("Fold", i, "\n")

      train_id <- unlist(train_clinical_IDs[[i]])
      test_id <- unlist(test_clinical_IDs[[i]])
      # test_id <- setdiff(rownames(clinical_data), train_id)

      # return(train_id)

      # Using the rownames of the selected IDs for x also for y
      # NA values are no longer used this way
      x_train <-
        clinical_data[as.character(train_id), ] %>% drop_na()

      y_train <- inc[as.character(rownames(x_train))]

      x_test <- clinical_data[as.character(test_id), ] %>% drop_na()

      y_test <- inc[rownames(x_test)]

      weights <-
        ifelse(y_train == "One", table(y_train)[[2]] / table(y_train)[[1]], 1)

      # Building the fit
      set.seed(993)
      fit <- caret::train(
        x = x_train,
        y = y_train,
        method = "glmnet",
        metric = "myLogLoss",
        tuneLength = 20,
        weights = weights,
        maximize = FALSE,
        trControl = my_control,
        importance = TRUE
      )

      # Building the prediction

      pred <-
        predict(fit, x_test, s = "lambda.min", type = "prob")$One
      roc <-
        roc(
          response = y_test,
          predictor = pred,
          levels = c("Zero", "One")
        )
      auc <- auc(roc)
      perf_list[[i]] <- auc
      var_imp <- varImp(fit)$importance
      vari <- var_imp$Overall
      names(vari) <- rownames(var_imp)
      vari <- vari[order(vari, decreasing = T)]

      var[[i]] <- vari

      # return(var)

    }

    set.seed(993)
    #Select best features

    var_fil <- lapply(var, function(x)
      names(x[x > 0]))

    # return(var_fil)
    var_rra <- aggregateRanks(var_fil)

    # return(var_rra)
    var_rra$adjP <- var_rra$Score * 100
    # return(var_rra)
    var_rra$adjP <- p.adjust(var_rra$adjP, "fdr")
    # var_rra$adjP <- 0.04
    # return(var_rra)
    var_sel <- var_rra %>% dplyr::filter(adjP < 0.05)

    var_sel <- as.numeric(unname(unlist(select(var_sel, "Name"))))

    # return(features_sel)

    saveRDS(dplyr::select(clinical_data, all_of(var_sel)),
            "clinical_selected.rds")


    return("Elasticnet clinical done!")
  }
