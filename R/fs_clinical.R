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

    test <- as.character(unlist(append(train_clinical_IDs[1], test_clinical_IDs[1])))

    # return(sort(test))

    mnsi <- dplyr::filter(phenotype_IDs, rownames(phenotype_IDs) %in% test)
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
      summaryFunction = mnLogLoss,
      sampling = NULL,
      allowParallel = T
    )

    # return(my_control)

    perf_list = list()
    var = list()


    # length(resamples)
    for (i in 1:2) {
      cat("Fold", i, "\n")

      train_id <- train_clinical_IDs[[i]]
      test_id <- test_clinical_IDs[[i]]

      x_train <- clinical_data[,!is.na(match(colnames(clinical_data), train_id))]
      y_train <- inc[as.character(colnames(x_train))]
      x_test <- clinical_data[, !is.na(match(colnames(clinical_data), test_id))]
      y_test <- inc[as.character(colnames(x_test))]



      return(table(y_train))
      # return(colnames(x_test))



      # Maybe not needed, because the selected IDs are checked through by the
      # selection of the training and testing IDs
      if (any(colnames(x_train) != names(y_train)) | any(colnames(x_test) != names(y_test))){
        message("Samples in train and test sets do not match!")
        next
      }

      return(y_train)

    }


  }
