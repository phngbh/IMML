#' Feature Selection Transcriptomics Data
#'
#' @description The function for the feature selection for transcriptomics data.
#'
#' @param trainIDs Set of training IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param testIDs Set of testing IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param dataIDs A table holding all the clinical IDs with the respective IDs
#' for each modality.
#' @param phenotypeIDs A table holding the clinical IDs with a variable indicating,
#' if the disease occurred or not.
#' @param transcriptomicsData A table holding the data for the transcriptomics.
#' @param geneAnotation An annotation file for the transcriptomics data.
#' @param lowestLevelPathways The lowest level pathways
#' @param seed The possibility to change the seed for the function.
#'
#' @return Returns the transcriptomics data table with the selected features.
#'
#' @author Ulrich Asemann

fs_transcriptomics <-
  function(trainIDs,
           testIDs,
           dataIDs,
           phenotypeIDs,
           transcriptomicsData,
           geneAnotation,
           lowestLevelPathways,
           seed = 123) {
    # Getting the ID sets from the transcriptomics
    trainTranscriptomicsIDs <-
      trainIDs$`Training Feature Selection IDs`$Transcriptomics
    testTranscriptomicsIDs <-
      testIDs$`Testing Feature Selection IDs`$Transcriptomics

    # Frame of phenotypes, with all used transcriptomic IDs and their inc3 value
    samples <- unlist(trainTranscriptomicsIDs) %>% unique()
    tranIDs <-
      dataIDs[match(samples, dataIDs$Clinical),] %>% dplyr::select("Transcriptomics") %>% unlist()

    info <-
      phenotypeIDs[as.character(samples), "inc3", drop = FALSE]
    rownames(info) <- tranIDs

    # Creating a model matrix
    info <- info %>% transmute(inc3 = as.character(inc3))

    modmatrix <- model.matrix( ~ 0 + ., data = info)

    # Creating a list for the gsea
    gsea = list(
      edge = list(),
      ilmn = list(),
      auc = list(),
      pathway = list()
    )

    # seed
    set.seed(seed)

    # Going through each set of IDs
    for (i in 1:length(trainTranscriptomicsIDs)) {
      cat("Iter ", i, "\n")
      # Selecting the dataIDs with one set of the dataPartitioning
      clinIDs <- unlist(trainTranscriptomicsIDs[[i]])

      tmpTranscriptomicsIDs <-
        dataIDs %>% filter(Clinical %in% clinIDs) %>%
        dplyr::select(Transcriptomics) %>% as.list() %>% unlist()

      # Getting the temporary data needed for calculations
      tmpData <-
        transcriptomicsData[,!is.na(match(colnames(transcriptomicsData), tmpTranscriptomicsIDs))]

      cat("...DE analysis\n")

      # Temporary model matrix
      tmpMod <- modmatrix[colnames(tmpData), , drop = FALSE]

      # Testing
      if (any(table(tmpMod[, 1]) < 2) |
          any(table(tmpMod[, 2]) < 2)) {
        cat("......One of the factor levels has less than 2 observations => Stop!")
        next
      }

      # lmFit
      fit <- lmFit(tmpData, tmpMod)

      # contrast
      contrast <-
        makeContrasts(inc31 - inc30, levels = colnames(coef(fit)))

      # contrasts.fit
      tmp <- contrasts.fit(fit, contrast)

      tmp <- eBayes(tmp)

      # topTable function
      topde <-
        topTable(tmp, sort.by = "P", n = Inf) %>% mutate(Probe = rownames(.)) %>%
        mutate(Name = geneAnotation$symbol[match(.$Probe, geneAnotation$Probe_Id)],
               EntrezID = geneAnotation$EntrezID[match(.$Probe, geneAnotation$Probe_Id)])

      cat("...GSEA\n")

      # Adding Additional variables (each iteration)
      topdeTmp <-
        topde %>% dplyr::select(EntrezID, P.Value) %>% na.omit()

      topdeP <- tapply(topdeTmp$P.Value, topdeTmp$EntrezID, min)

      # Initializing the ranklist
      ranklist <- vector(mode = "numeric", length = length(topdeP))
      names(ranklist) = names(topdeP)

      # Temporary data frame with all values
      tmp <-
        filter(topde, EntrezID %in% names(topdeP) &
                 P.Value %in% topdeP)

      # Adding values to the ranklist vector
      tTmp <- tmp %>% dplyr::select(EntrezID, t)
      ranklist[as.character(tTmp$EntrezID)] <- tTmp$t

      # Initializing the genelistTmp variable
      genelistTmp <-
        tmp %>% dplyr::select(Probe, EntrezID, t, Name)
      rownames(genelistTmp) <- as.character(1:nrow(genelistTmp))

      # Sorting the ranklist
      ranklist <- sort(ranklist)


      # Muting warning messages for now
      suppressMessages({
        suppressWarnings({
          # genesetReactome
          genesetReactome <- reactomePathways(names(ranklist))
          genesetReactome <-
            genesetReactome[intersect(names(genesetReactome), lowestLevelPathways)]

          # Set seed
          set.seed(seed)

          # fgsea function
          fgseaTmp <- fgsea(
            pathways = genesetReactome,
            stats = ranklist,
            minSize = 15,
            maxSize = 200,
            eps = 0
          ) %>% arrange(pval) %>% filter(padj < 0.1)
        })
      })

      # Test
      if (nrow(fgseaTmp) == 0) {
        print("No significant pathway found")
        next
      }

      # Defining edgeTmp
      edgeTmp <- fgseaTmp$leadingEdge %>% unlist() %>% unique()

      # Add values to list
      gsea$edge[[i]] <- edgeTmp
      gsea$ilmn[[i]] <-
        genelistTmp$Probe[match(edgeTmp, genelistTmp$EntrezID)]
      gsea$pathway[[i]] <- fgseaTmp

      cat("...Elastic net\n")

      # Test
      if (is.null(gsea$ilmn[[i]])) {
        print("No significant pathway found")
        next
      }

      # Transcriptomic IDs used for training
      resamples <-
        transcriptomicsData[, colnames(transcriptomicsData) %in% tmpTranscriptomicsIDs] %>%
        colnames()

      # Training/Testing
      probelistTmp <- gsea$ilmn[[i]]
      testID <-
        colnames(transcriptomicsData[,!(colnames(transcriptomicsData) %in% resamples)])
      xTrain <-
        transcriptomicsData[probelistTmp, resamples] %>% t()
      yTrain <- info[resamples, "inc3"]
      yTrain <-
        ifelse(yTrain == 1, "One", "Zero") %>% factor(levels = c("One", "Zero"))
      xTest <- transcriptomicsData[probelistTmp, testID] %>% t()
      yTest <- info[testID, "inc3"]
      yTest <-
        ifelse(yTest == 1, "One", "Zero") %>% factor(levels = c("One", "Zero"))

      # Control
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

      # Set seed
      set.seed(seed)

      # Create fit
      fit <- train(
        x = xTrain,
        y = yTrain,
        method = "glmnet",
        metric = "ROC",
        tuneLength = 20,
        maximize = TRUE,
        trControl = myControl,
        importance = TRUE
      )

      # Prediction
      pred <- predict(fit,
                      xTest,
                      s = "lambda.min",
                      type = "prob")$One

      roc <- roc(
        response = yTest,
        predictor = pred,
        levels = c("Zero", "One")
      )

      auc <- auc(roc)

      # Adding to list
      gsea$auc[[i]] <- auc
    }

    # saveRDS(gsea, "gsea_test_list.rds")

    # Select the top significant pathways
    keep <- vector("logical", length = length(gsea$auc))

    for (i in 1:length(keep)) {
      if (is.null(gsea$auc[[i]])) {
        next
      } else if (gsea$auc[[i]] <= 0.5) {
        next
      } else {
        keep[i] <- TRUE
      }
    }

    pwlist <- gsea$pathway[keep]

    pathways <- list(pw <- list(),
                     pval <- list(),
                     NES <- list())
    for (i in 1:length(pwlist)) {
      if (is.null(pwlist[[i]])) {
        next
      }
      pathways$pw[[i]] <- pwlist[[i]]$pathway
      pathways$pval[[i]] <- pwlist[[i]]$pval
      pathways$NES[[i]] <- pwlist[[i]]$NES
    }

    pwlistAgg <- aggregateRanks(pathways$pw)
    pwlistAgg$adjP <- pwlistAgg$Score * length(pathways$pw)
    pwlistAgg$adjP <- p.adjust(pwlistAgg$adjP, method = "fdr")

    toppw <- rownames(filter(pwlistAgg, adjP < 0.05))

    # Final gene set enrichment analysis on the selected pathways
    fit <-
      lmFit(transcriptomicsData[,!is.na(match(colnames(transcriptomicsData),
                                              rownames(modmatrix))), drop = F],
            modmatrix[!is.na(match(rownames(modmatrix), colnames(transcriptomicsData))), , drop = F])

    contrast <-
      makeContrasts(inc31 - inc30, levels = colnames(coef(fit)))


    tmp <- contrasts.fit(fit, contrast)
    tmp <- eBayes(tmp)

    topde <- topTable(tmp,
                      sort.by = "P",
                      n = Inf)
    topde <- topde %>%
      mutate(Gene = rownames(topde)) %>%
      mutate(Name = geneAnotation$symbol[match(.$Gene, geneAnotation$Probe_Id)],
             EntrezID = geneAnotation$EntrezID[match(.$Gene, geneAnotation$Probe_Id)])

    topdeTmp <-
      topde %>% dplyr::select(EntrezID, P.Value) %>% na.omit()

    topdeP <- tapply(topdeTmp$P.Value, topdeTmp$EntrezID, min)

    ranklist <- vector(mode = "numeric",
                       length = length(topdeP))
    names(ranklist) <- names(topdeP)

    # Temporary data frame with all values
    tmp <-
      filter(topde, EntrezID %in% names(topdeP) &
               P.Value %in% topdeP)

    # Adding values to the ranklist vector
    tTmp <- tmp %>% dplyr::select(EntrezID, t)
    ranklist[as.character(tTmp$EntrezID)] <- tTmp$t

    # Initializing the genelist variable
    genelist <- tmp %>% dplyr::select(Gene, EntrezID, t, Name)
    rownames(genelist) <- as.character(1:nrow(genelist))

    # Sorting the ranklist
    ranklist <- sort(ranklist)

    # Geneset reactome
    # Muteing warning messages for now
    suppressMessages({
      suppressWarnings({
        # genesetReactome
        genesetReactome <- reactomePathways(names(ranklist))
        genesetReactome <-
          genesetReactome[intersect(names(genesetReactome), toppw)]

        # Set seed
        set.seed(seed)

        # fgsea function
        fgseaRes <- fgsea(
          pathways = genesetReactome,
          stats = ranklist,
          minSize = 15,
          maxSize = 200
        ) %>% arrange(pval) #%>% filter(padj < 0.1)

        # return(fgseaRes)

        # Adding values to leadingEdge
        fgseaRes <-
          fgseaRes %>% mutate(
            leadingEdge = mapIdsList(
              x = org.Hs.eg.db,
              keys = leadingEdge,
              keytype = "ENTREZID",
              column = "SYMBOL"
            )
          )

        # Saving the result to the enviroment
        saveRDS(fgseaRes, "trans_gsea_final.rds")

        # Extract selected features for training
        edge <- fgseaRes$leadingEdge %>% unlist() %>% unique()
        edgeEntrez <-
          AnnotationDbi::select(
            org.Hs.eg.db,
            keys = edge,
            columns = "ENTREZID",
            keytype = "SYMBOL"
          )
        probe <-
          genelist$Gene[genelist$EntrezID %in% edgeEntrez$ENTREZID]
        dataSelected <- transcriptomicsData[probe, , drop = FALSE]

        # Saving the data in a new frame
        saveRDS(dataSelected, "transcriptomics_selected.rds")

      })
    })

    return("Feature selection transcriptomics done!")
  }
