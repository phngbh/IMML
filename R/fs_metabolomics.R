#' Feature Selection for Metabolomics
#'
#' @param train_IDs
#' @param test_IDs
#' @param data_IDs
#' @param phenotype_IDs
#' @param metabolomics_data
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
fs_metabolomics <-
  function(train_IDs,
           test_IDs,
           data_IDs,
           phenotype_IDs,
           metabolomics_data,
           pathway_list,
           low_lev_path,
           anno_fil,
           seed = 123) {
    # Selecting the IDs
    train_metabolomics_IDs <-
      train_IDs$`Feature Selection IDs`$Metabolomics
    test_metabolomics_IDs <-
      test_IDs$`Feature Selection IDs`$Metabolomics

    # Data frame of phenotypes with all used IDs and their inc3 value
    samples <- unlist(train_metabolomics_IDs) %>% unique()
    meta_IDs <-
      data_IDs[match(samples, data_IDs$Clinical), ] %>% dplyr::select("Metabolomics") %>% unlist()

    info <-
      phenotype_IDs[as.character(samples), "inc3", drop = FALSE]
    rownames(info) <- meta_IDs

    # Creating a model matrix
    info <- info %>% transmute(inc3 = as.character(inc3))

    modmatrix <- model.matrix(~ 0 + ., data = info)

    # Creating a list for the gsea
    gsea = list(
      edge = list(),
      probe = list(),
      auc = list(),
      pathway = list()
    )

    # Setting the seed for the methode
    set.seed(seed)

    # Going through each set of IDs

    # normal length of iteration
    # i in 1:length(train_metabolomics_IDs)

    # Adding parallel computation
    # Removed for now
    # cl <- makeCluster(2)
    # registerDoSNOW(cl)
    #
    # foreach (i = 1:100, .packages = c("org.Hs.eg.db",
    #                                   "dplyr",
    #                                   "tidyr",
    #                                   "caret",
    #                                   "BiocManager",
    #                                   "limma",
    #                                   "ChAMP",
    #                                   "fgsea",
    #                                   "glmnet",
    #                                   "themis",
    #                                   "pROC",
    #                                   "RobustRankAggreg",
    #                                   "knitr")) %dopar%
    for (i in 1:100) {
      cat("Iter ", i, "\n")

      # Selecting the data_IDs with one set of the data_partitioning
      clin_IDs <- unlist(train_metabolomics_IDs[[i]])

      tmp_metabolomics_IDs <-
        data_IDs %>% filter(Clinical %in% clin_IDs) %>%
        dplyr::select(Metabolomics) %>% as.list() %>% unlist()

      # Getting the temporary data needed for calculations
      tmp_data <-
        metabolomics_data[, !is.na(match(colnames(metabolomics_data), tmp_metabolomics_IDs))]

      cat("...DE analysis\n")

      # Temporary model matrix
      tmp_mod <- modmatrix[colnames(tmp_data), , drop = FALSE]

      # Testing
      if (any(table(tmp_mod[, 1]) < 2) |
          any(table(tmp_mod[, 2]) < 2)) {
        cat("......One of the factor levels has less than 2 observations => Stop!")
        next
      }

      # lmFit
      fit <- lmFit(tmp_data, tmp_mod)

      # contrast
      contrast <-
        makeContrasts(inc31 - inc30, levels = colnames(coef(fit)))

      # contrasts.fit
      tmp <- contrasts.fit(fit, contrast)

      tmp <- eBayes(tmp)

      #topTable function
      topde <-
        topTable(tmp, sort.by = "P", n = Inf) %>% mutate(ID = rownames(.)) %>%
        mutate(Name = anno_fil$Biochemical_F4[match(.$ID, anno_fil$Mnumber)],
               ChEBI = anno_fil$ChEBI[match(.$ID, anno_fil$Mnumber)])

      cat("...GSEA\n")

      # ranklist
      ranklist <- topde$t
      names(ranklist) <- topde$ChEBI
      ranklist <- sort(ranklist)

      # set seed
      set.seed(seed)

      # gsea
      suppressMessages({
        suppressWarnings({
          fgsea_tmp <- fgsea(
            pathways = pathway_list,
            stats = ranklist,
            minSize = 1,
            maxSize = 200,
            eps = 0
          ) %>% dplyr::arrange(pval) %>% dplyr::filter(padj < 0.1)
        })
      })
      # Test
      if (nrow(fgsea_tmp) == 0) {
        message("No significant pathway found")
        next
      }

      # Defining edge_tmp
      edge_tmp <- fgsea_tmp$leadingEdge %>% unlist() %>% unique()

      # Adding values to list
      gsea$edge[[i]] <- edge_tmp
      gsea$probe[[i]] <-
        anno_fil$Mnumber[match(edge_tmp, anno_fil$ChEBI)]
      gsea$pathway[[i]] <- fgsea_tmp

      # return(fgsea_tmp)

      cat("...Elastic net\n")

      # Metabolomic IDs used for training
      resamples <-
        metabolomics_data[, colnames(metabolomics_data) %in% tmp_metabolomics_IDs] %>%
        colnames()

      # Training/Testing
      probelist_tmp <- gsea$probe[[i]]
      test_ID <-
        colnames(metabolomics_data[, !(colnames(metabolomics_data) %in% resamples)])
      x_train <- metabolomics_data[probelist_tmp, resamples] %>% t()
      y_train <- info[resamples, "inc3"]
      y_train <-
        ifelse(y_train == 1, "One", "Zero") %>% factor(levels = c("One", "Zero"))
      x_test <- metabolomics_data[probelist_tmp, test_ID] %>% t()
      y_test <- info[test_ID, "inc3"]
      y_test <-
        ifelse(y_test == 1, "One", "Zero") %>% factor(levels = c("One", "Zero"))

      # My control
      my_control <- trainControl(
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
        x = x_train,
        y = y_train,
        method = "glmnet",
        metric = "ROC",
        tuneLength = 20,
        maximize = T,
        trControl = my_control,
        importance = TRUE
      )

      # Predict
      pred = predict(fit, x_test, s = "lambda.min", type = "prob")$One
      roc <-
        roc(
          response = y_test,
          predictor = pred,
          levels = c("Zero", "One")
        )
      auc = auc(roc)
      gsea$auc[[i]] = auc

    }

    # Stop the cluster
    # stopCluster(cl)

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

    # return(keep)

    # Pathways
    pwlist <- gsea$pathway[keep]

    # return(pwlist)

    pathways <- list(pw = list(),
                     pval = list(),
                     NES = list())

    for (i in 1:length(pwlist)) {
      if (is.null(pwlist[[i]])) {
        next
      }
      pathways$pw[[i]] = pwlist[[i]]$pathway
      pathways$pval[[i]] = pwlist[[i]]$pval
      pathways$NES[[i]] = pwlist[[i]]$NES
    }

    pwlist_agg <- aggregateRanks(pathways$pw)

    # return(pwlist_agg)

    pwlist_agg$adjP <- pwlist_agg$Score * length(pathways$pw)

    # return(pwlist_agg)
    # return(length(pathways$pw))

    pwlist_agg$adjP <- p.adjust(pwlist_agg$adjP, method = "fdr")

    # return(pwlist_agg)

    toppw <- rownames(filter(pwlist_agg, adjP < 0.05))

    # is empty because of wrong p-values
    # return(toppw)

    # Final gene set enrichment analysis on the selected pathways
    # Fit
    fit <-
      lmFit(metabolomics_data[,!is.na(match(colnames(metabolomics_data),
                                            rownames(modmatrix))), drop = F],
            modmatrix[!is.na(match(rownames(modmatrix), colnames(metabolomics_data))), , drop = F])

    # Contrast
    contrast <-
      makeContrasts(inc31 - inc30, levels = colnames(coef(fit)))
    tmp <- contrasts.fit(fit, contrast)
    tmp <- eBayes(tmp)

    # topTable
    topde <-
      topTable(tmp, sort.by = "P", n = Inf) %>% mutate(ID = rownames(.)) %>%
      mutate(Name = anno_fil$Biochemical_F4[match(.$ID, anno_fil$Mnumber)],
             ChEBI = anno_fil$ChEBI[match(.$ID, anno_fil$Mnumber)])

    # return(topde)

    # ranklist
    ranklist <- topde$t
    topde$ChEBI <- anno_fil$ChEBI[match(topde$ID, anno_fil$Mnumber)]
    names(ranklist) <- topde$ChEBI
    ranklist <- sort(ranklist)

    suppressMessages({
      suppressWarnings({
        # Pathways
        geneset_reactome <- reactomePathways(names(ranklist))
        geneset_reactome <-
          geneset_reactome[intersect(names(geneset_reactome), toppw)]


        # Seed
        set.seed(seed)


        # Change pathway, either pathway_list or geneset_reactome
        # Problem with minSize, only returns something when size is 1

        # fgsea
        fgseaRes <- fgsea(
          pathways = geneset_reactome,
          stats    = ranklist,
          minSize  = 1,
          maxSize  = 200
        ) %>% arrange(pval) #%>% filter(padj < 0.3)

        # return(fgseaRes)

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

        # Extract selected features for training
        edge <- fgseaRes$leadingEdge %>% unlist() %>% unique()
        probe <- anno_fil$Mnumber[anno_fil$ChEBI %in% edge]
        dat_selected <- metabolomics_data[probe, , drop = FALSE]

        saveRDS(dat_selected, "metabolomics_selected.rds")

      })
    })

    # End of function
    return("Feature selection metabolomics done!")

  }
