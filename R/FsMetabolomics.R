#' Feature Selection Metabolomics Data
#'
#' @description The function for the feature selection for metabolomics data.
#'
#' @param trainIDs Set of training IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param testIDs Set of testing IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param dataIDs A table holding all the clinical IDs with the respective IDs
#' for each modality.
#' @param phenotypeIDs A table holding the clinical IDs with a variable indicating,
#' if the disease occurred or not.
#' @param metabolomicsData A table holding the data for the metabolomics.
#' @param seed The possibility to change the seed for the function.
#' @param pathwayList A list containing the pathways for the metabolomics.
#' @param lowestLevelPathways The lowest level pathways.
#' @param annotationFile An annotation file containing additional information for the metabolomics.
#'
#' @return Returns the metabolomics data table with the selected features.
#'
#' @author Ulrich Asemann

FsMetabolomics <-
  function(trainIDs,
           testIDs,
           dataIDs,
           phenotypeIDs,
           metabolomicsData,
           pathwayList,
           lowestLevelPathways,
           annotationFile,
           seed = 123) {
    # Selecting the IDs
    trainMetabolomicsIDs <-
      trainIDs$`Training Feature Selection IDs`$Metabolomics
    testMetabolomicsIDs <-
      testIDs$`Testing Feature Selection IDs`$Metabolomics

    # Data frame of phenotypes with all used IDs and their inc3 value
    samples <- unlist(trainMetabolomicsIDs) %>% unique()
    metabolomicsIDs <-
      dataIDs[match(samples, dataIDs$Clinical),] %>% dplyr::select("Metabolomics") %>% unlist()

    info <-
      phenotypeIDs[as.character(samples), "inc3", drop = FALSE]
    rownames(info) <- metabolomicsIDs

    # Creating a model matrix
    info <- info %>% transmute(inc3 = as.character(inc3))

    modmatrix <- model.matrix( ~ 0 + ., data = info)

    # Creating a list for the gsea
    gsea = list(
      edge = list(),
      probe = list(),
      auc = list(),
      pathway = list()
    )

    # Setting the seed for the function
    set.seed(seed)

    # Going through each set of IDs
    # 1:length(trainMetabolomicsIDs)
    for (i in 1:length(trainMetabolomicsIDs)) {
      cat("Iter ", i, "\n")

      # Selecting the dataIDs from one set of the dataPartitioning()
      clinicalIDs <- unlist(trainMetabolomicsIDs[[i]])

      tmpMetabolomicsIDs <-
        dataIDs %>% filter(Clinical %in% clinicalIDs) %>%
        dplyr::select(Metabolomics) %>% as.list() %>% unlist()

      # Getting the temporary data needed for calculations
      tmpData <-
        metabolomicsData[,!is.na(match(colnames(metabolomicsData), tmpMetabolomicsIDs))]

      cat("...DE analysis\n")

      # Temporary model matrix
      tmpModmatrix <- modmatrix[colnames(tmpData), , drop = FALSE]

      # Testing
      if (any(table(tmpModmatrix[, 1]) < 2) |
          any(table(tmpModmatrix[, 2]) < 2)) {
        cat("......One of the factor levels has less than 2 observations => Stop!")
        next
      }

      # lmFit
      fit <- lmFit(tmpData, tmpModmatrix)

      # contrast
      contrast <-
        makeContrasts(inc31 - inc30, levels = colnames(coef(fit)))

      # contrasts.fit
      tmp <- contrasts.fit(fit, contrast)

      tmp <- eBayes(tmp)

      # topTable function
      topde <-
        topTable(tmp, sort.by = "P", n = Inf) %>% mutate(ID = rownames(.)) %>%
        mutate(Name = annotationFile$Biochemical_F4[match(.$ID, annotationFile$Mnumber)],
               ChEBI = annotationFile$ChEBI[match(.$ID, annotationFile$Mnumber)])

      cat("...GSEA\n")

      # ranklist
      ranklist <- topde$t
      names(ranklist) <- topde$ChEBI
      ranklist <- sort(ranklist)

      # set seed
      set.seed(seed)

      # Fast gsea
      suppressMessages({
        suppressWarnings({
          fgseaTmp <- fgsea(
            pathways = pathwayList,
            stats = ranklist,
            minSize = 1,
            maxSize = 200,
            eps = 0
          ) %>% dplyr::arrange(pval) %>% dplyr::filter(padj < 0.1)
        })
      })

      # Test
      if (nrow(fgseaTmp) == 0) {
        message("No significant pathway found")
        next
      }

      # Defining edgeTmp
      edgeTmp <- fgseaTmp$leadingEdge %>% unlist() %>% unique()

      # Adding values to list
      gsea$edge[[i]] <- edgeTmp
      gsea$probe[[i]] <-
        annotationFile$Mnumber[match(edgeTmp, annotationFile$ChEBI)]
      gsea$pathway[[i]] <- fgseaTmp

      cat("...Elastic net\n")

      # Metabolomic IDs used for training
      resamples <-
        metabolomicsData[, colnames(metabolomicsData) %in% tmpMetabolomicsIDs] %>%
        colnames()

      # Training/Testing
      probelistTmp <- gsea$probe[[i]]
      testID <-
        colnames(metabolomicsData[,!(colnames(metabolomicsData) %in% resamples)])
      xTrain <- metabolomicsData[probelistTmp, resamples] %>% t()
      yTrain <- info[resamples, "inc3"]
      yTrain <-
        ifelse(yTrain == 1, "One", "Zero") %>% factor(levels = c("One", "Zero"))
      xTest <- metabolomicsData[probelistTmp, testID] %>% t()
      yTest <- info[testID, "inc3"]
      yTest <-
        ifelse(yTest == 1, "One", "Zero") %>% factor(levels = c("One", "Zero"))

      # My control
      myControl <- trainControl(
        method = "repeatedcv",
        number = 5,
        repeats = 2,
        savePredictions = "final",
        classProbs = TRUE,
        summaryFunction = twoClassSummary,
        sampling = "smote",
        allowParallel = TRUE
      )

      # Seed
      set.seed(seed)

      # Fit
      fit <- caret::train(
        x = xTrain,
        y = yTrain,
        method = "glmnet",
        metric = "ROC",
        tuneLength = 20,
        maximize = T,
        trControl = myControl,
        importance = TRUE
      )

      # Predict
      pred = predict(fit, xTest, s = "lambda.min", type = "prob")$One
      roc <-
        roc(
          response = yTest,
          predictor = pred,
          levels = c("Zero", "One")
        )
      auc = auc(roc)
      gsea$auc[[i]] = auc

    }

    # Save results
    saveRDS(gsea, "meta_gsea_list.rds")

    # Select the top significant pathways
    keep <- vector("logical", length = length(gsea$auc))
    for (i in 1:length(keep)) {
      if (is.null(gsea$auc[[i]])) {
        next
      } else if (gsea$auc[[i]] <= 0.5) {
        next
      } else{
        keep[i] <- TRUE
      }
    }

    # Pathways
    pathwaysList <- gsea$pathway[keep]

    pathways <- list(pw = list(),
                     pval = list(),
                     NES = list())

    for (i in 1:length(pathwaysList)) {
      if (is.null(pathwaysList[[i]])) {
        next
      }
      pathways$pw[[i]] = pathwaysList[[i]]$pathway
      pathways$pval[[i]] = pathwaysList[[i]]$pval
      pathways$NES[[i]] = pathwaysList[[i]]$NES
    }

    pathwaysListAggregated <- aggregateRanks(pathways$pw)

    pathwaysListAggregated$adjP <-
      pathwaysListAggregated$Score * length(pathways$pw)

    pathwaysListAggregated$adjP <-
      p.adjust(pathwaysListAggregated$adjP, method = "fdr")

    toppw <- rownames(filter(pathwaysListAggregated, adjP < 0.05))

    # Final gene set enrichment analysis on the selected pathways
    # Fit
    fit <-
      lmFit(metabolomicsData[, !is.na(match(colnames(metabolomicsData),
                                            rownames(modmatrix))), drop = F],
            modmatrix[!is.na(match(rownames(modmatrix), colnames(metabolomicsData))), , drop = F])

    # Contrast
    contrast <-
      makeContrasts(inc31 - inc30, levels = colnames(coef(fit)))
    tmp <- contrasts.fit(fit, contrast)
    tmp <- eBayes(tmp)

    # topTable
    topde <-
      topTable(tmp, sort.by = "P", n = Inf) %>% mutate(ID = rownames(.)) %>%
      mutate(Name = annotationFile$Biochemical_F4[match(.$ID, annotationFile$Mnumber)],
             ChEBI = annotationFile$ChEBI[match(.$ID, annotationFile$Mnumber)])

    # ranklist
    ranklist <- topde$t
    topde$ChEBI <-
      annotationFile$ChEBI[match(topde$ID, annotationFile$Mnumber)]
    names(ranklist) <- topde$ChEBI
    ranklist <- sort(ranklist)

    suppressMessages({
      suppressWarnings({
        # Pathways
        genesetReactome <- reactomePathways(names(ranklist))
        genesetReactome <-
          genesetReactome[intersect(names(genesetReactome), toppw)]

        # Seed
        set.seed(seed)

        # fgsea
        fgseaRes <- fgsea(
          pathways = genesetReactome,
          stats    = ranklist,
          # minSize  = 1,
          minSize = 5,
          maxSize  = 200
        ) %>% arrange(pval) #%>% filter(padj < 0.3)

        fgseaRes <-
          fgseaRes %>% mutate(
            leadingEdge = mapIdsList(
              x = org.Hs.eg.db,
              keys = leadingEdge,
              keytype = "ENTREZID",
              column = "SYMBOL"
            )
          )

        saveRDS(fgseaRes, "meta_gsea_final.rds")

        # Extract selected features
        edge <- fgseaRes$leadingEdge %>% unlist() %>% unique()
        probe <-
          annotationFile$Mnumber[annotationFile$ChEBI %in% edge]
        dataSelected <- metabolomicsData[probe, , drop = FALSE]

        # Save results
        saveRDS(dataSelected, "metabolomics_selected.rds")

      })
    })

    # End of function
    return("Feature selection metabolomics done!")

  }
