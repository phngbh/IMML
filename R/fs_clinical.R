#' Title
#'
#' @param train_IDs
#' @param test_IDs
#' @param data_IDs
#'
#' @return
#' @export
#'
#' @examples
fs_clinical <-
  function(train_IDs,
           test_IDs,
           data_IDs,
           phenotype_IDs,
           clinical_data) {
    # Getting IDs
    train_clinical_IDs <-
      train_IDs$`Feature Selection IDs`$Clinical
    test_clinical_IDs <-
      test_IDs$`Feature Selection IDs`$Clinical

    #

    # # return(clinical_data[,as.numeric("2")])
    #
    # test <-
    #   as.character(unlist(append(
    #     train_clinical_IDs[1], test_clinical_IDs[1]
    #   )))
    #
    # # return(sort(test))
    #
    # mnsi <-
    #   dplyr::filter(phenotype_IDs, rownames(phenotype_IDs) %in% test)
    # inc <- mnsi$inc3
    # names(inc) <- rownames(mnsi)
    # inc <-
    #   ifelse(inc == 1, "One", "Zero") %>% factor(levels = c("One", "Zero"))
    #
    # # return(length(inc))

    mnsi <- filter(phenotype_IDs,!is.na(inc3))

    inc <- mnsi$inc3

    names(inc) <- rownames(mnsi)

    inc <-
      ifelse(inc == 1, "One", "Zero") %>% factor(levels = c("One", "Zero"))

    # return(inc)






    #Do elastic net
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
    for (i in 1:20) {
      cat("Fold", i, "\n")

      train_id <- unlist(train_clinical_IDs[[i]])
      test_id <- unlist(test_clinical_IDs[[i]])
      # test_id <- setdiff(rownames(clinical_data), train_id)

      # return(test_id)

      # Workaround at the moment, because of missing data in clinical_processed
      x_train <-
        clinical_data[as.character(train_id), ] %>% drop_na()
      y_train <- inc[as.character(train_id)]
      x_test <- clinical_data[as.character(test_id), ] #%>% drop_na()
      y_test <- inc[rownames(x_test)]

      return(rownames(x_train))

      weights <-
        ifelse(y_train == "One", table(y_train)[[2]] / table(y_train)[[1]], 1)

      # return(weights)

      # return(any(rownames(x_train) != names(y_train)))
      # return(rownames(x_train))
      # return(names(y_train))
      # return(colnames(x_test))

      # Maybe not needed, because the selected IDs are checked through by the
      # selection of the training and testing IDs
      if (any(rownames(x_train) != names(y_train)) |
          any(rownames(x_test) != names(y_test))) {
        message("Samples in train and test sets do not match!")
        next
      }

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

      # return((fit))

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

      # return(var_fil)
      #var = var[var != 0]
      var[[i]] <- vari

      # return(vari)

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

    # return(rownames(var_sel))

    return(clinical_data[, as.numeric(rownames(var_sel)) %in% 1:ncol(clinical_data)])

    saveRDS(clinical_data[, as.numeric(rownames(var_sel))], "clinical_selected.rds")


    return("done")
  }
