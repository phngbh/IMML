#' Feature Selection for Transcriptomics
#'
#' @param data_IDs
#' @param transcriptomics_data
#' @param train_IDs
#' @param test_IDs
#' @param phenotype_IDs
#' @param gene_anotation
#' @param low_lev_path
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
fs_transcriptomics <-
  function(train_IDs,
           test_IDs,
           data_IDs,
           phenotype_IDs,
           transcriptomics_data,
           gene_anotation,
           low_lev_path,
           seed = 123) {
    # Working with gene_IDs, file from Phong, with annotation for the transcriptomics
    # Ignored annotation_human_..., maybe implement later
    #
    # Getting the ID sets from the transcriptomics
    train_transcriptomics_IDs <-
      train_IDs$`Feature Selection IDs`$Transcriptomics
    test_transcriptomics_IDs <-
      test_IDs$`Feature Selection IDs`$Transcriptomics

    # Frame of phenotypes, with all used transcriptomic IDs and their inc3 value
    samples <- unlist(train_transcriptomics_IDs) %>% unique()
    tran_IDs <-
      data_IDs[match(samples, data_IDs$Clinical),] %>% dplyr::select("Transcriptomics") %>% unlist()

    info <-
      phenotype_IDs[as.character(samples), "inc3", drop = FALSE]
    rownames(info) <- tran_IDs


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

    # Setting the seed for the methode
    set.seed(seed)

    # Going through each set of IDs

    # Normal length of Iteration:
    # i in 1:length(train_transcriptomics_IDs)


    for (i in 1:length(train_transcriptomics_IDs)) {
      cat("Iter ", i, "\n")
      # Selecting the data_IDs with one set of the data_partitioning
      clin_IDs <- unlist(train_transcriptomics_IDs[[i]])

      tmp_transcriptomics_IDs <-
        data_IDs %>% filter(Clinical %in% clin_IDs) %>%
        dplyr::select(Transcriptomics) %>% as.list() %>% unlist()

      # Getting the temporary data needed for calculations
      tmp_data <-
        transcriptomics_data[,!is.na(match(colnames(transcriptomics_data), tmp_transcriptomics_IDs))]

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

      # topTable function
      topde <-
        topTable(tmp, sort.by = "P", n = Inf) %>% mutate(Probe = rownames(.)) %>%
        mutate(Name = gene_anotation$symbol[match(.$Probe, gene_anotation$Probe_Id)],
               EntrezID = gene_anotation$EntrezID[match(.$Probe, gene_anotation$Probe_Id)])

      cat("...GSEA\n")

      # return(summary(topde))

      # Adding Additional variables (each iteration)
      topde_tmp <-
        topde %>% dplyr::select(EntrezID, P.Value) %>% na.omit()

      topde_p <- tapply(topde_tmp$P.Value, topde_tmp$EntrezID, min)

      # Initializing the ranklist
      ranklist <- vector(mode = "numeric", length = length(topde_p))
      names(ranklist) = names(topde_p)

      # Temporary data frame with all values
      tmp <-
        filter(topde, EntrezID %in% names(topde_p) &
                 P.Value %in% topde_p)

      # Adding values to the ranklist vector
      t_tmp <- tmp %>% dplyr::select(EntrezID, t)
      ranklist[as.character(t_tmp$EntrezID)] <- t_tmp$t

      # Initializing the genelist_tmp variable
      genelist_tmp <-
        tmp %>% dplyr::select(Probe, EntrezID, t, Name)
      rownames(genelist_tmp) <- as.character(1:nrow(genelist_tmp))

      # Sorting the ranklist
      ranklist <- sort(ranklist)


      # Muteing warning messages for now
      suppressMessages({
        suppressWarnings({
          # geneset_reactome
          geneset_reactome <- reactomePathways(names(ranklist))
          geneset_reactome <-
            geneset_reactome[intersect(names(geneset_reactome), low_lev_path)]

          # Set seed
          set.seed(seed)

          # fgsea function
          fgsea_tmp <- fgsea(
            pathways = geneset_reactome,
            stats = ranklist,
            minSize = 15,
            maxSize = 200,
            eps = 0
          ) %>% arrange(pval) %>% filter(padj < 0.1)
        })
      })

      # Test
      if (nrow(fgsea_tmp) == 0) {
        print("No significant pathway found")
        next
      }

      # Defining edge_tmp
      edge_tmp <- fgsea_tmp$leadingEdge %>% unlist() %>% unique()

      # Add values to list
      gsea$edge[[i]] <- edge_tmp
      gsea$ilmn[[i]] <-
        genelist_tmp$Probe[match(edge_tmp, genelist_tmp$EntrezID)]
      gsea$pathway[[i]] <- fgsea_tmp

      cat("...Elastic net\n")

      # Test
      if (is.null(gsea$ilmn[[i]])) {
        print("No significant pathway found")
        next
      }

      # Transcriptomic IDs used for training
      resamples <-
        transcriptomics_data[, colnames(transcriptomics_data) %in% tmp_transcriptomics_IDs] %>%
        colnames()

      # Training/Testing
      probelist_tmp <- gsea$ilmn[[i]]
      test_ID <-
        colnames(transcriptomics_data[,!(colnames(transcriptomics_data) %in% resamples)])
      x_train <-
        transcriptomics_data[probelist_tmp, resamples] %>% t()
      y_train <- info[resamples, "inc3"]
      y_train <-
        ifelse(y_train == 1, "One", "Zero") %>% factor(levels = c("One", "Zero"))
      x_test <- transcriptomics_data[probelist_tmp, test_ID] %>% t()
      y_test <- info[test_ID, "inc3"]
      y_test <-
        ifelse(y_test == 1, "One", "Zero") %>% factor(levels = c("One", "Zero"))

      # Control
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

      # Set seed
      set.seed(seed)

      # Create fit
      fit <- train(
        x = x_train,
        y = y_train,
        method = "glmnet",
        metric = "ROC",
        tuneLength = 20,
        maximize = TRUE,
        trControl = my_control,
        importance = TRUE
      )

      # Prediction
      pred <- predict(fit,
                      x_test,
                      s = "lambda.min",
                      type = "prob")$One

      roc <- roc(
        response = y_test,
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

    pwlist_agg <- aggregateRanks(pathways$pw)
    pwlist_agg$adjP <- pwlist_agg$Score * length(pathways$pw)
    pwlist_agg$adjP <- p.adjust(pwlist_agg$adjP, method = "fdr")

    toppw <- rownames(filter(pwlist_agg, adjP < 0.05))


    # # Dead Code!!!
    #
    # toppw_pval <- list()
    # for (p in toppw) {
    #   tmplist <- list()
    #
    #   for (i in 1:length(pwlist)) {
    #     ind <- which(pathways$pw[[i]] == p)
    #
    #     if (length(ind) > 0) {
    #       tmplist[[i]] <- pathways$pval[[i]][ind]
    #     } else {
    #       tmplist[[i]] <- NULL
    #     }
    #
    #   }
    #   tmpvec <- unlist(tmplist)
    #   toppw_pval[[p]] <- mean(tmpvec)
    # }

    # Final gene set enrichment analysis on the selected pathway
    fit <-
      lmFit(transcriptomics_data[,!is.na(match(colnames(transcriptomics_data),
                                               rownames(modmatrix))), drop = F],
            modmatrix[!is.na(match(rownames(modmatrix), colnames(transcriptomics_data))), , drop = F])

    contrast <-
      makeContrasts(inc31 - inc30, levels = colnames(coef(fit)))


    tmp <- contrasts.fit(fit, contrast)
    tmp <- eBayes(tmp)

    topde <- topTable(tmp,
                      sort.by = "P",
                      n = Inf)
    topde <- topde %>%
      mutate(Gene = rownames(topde)) %>%
      mutate(Name = gene_anotation$symbol[match(.$Gene, gene_anotation$Probe_Id)],
             EntrezID = gene_anotation$EntrezID[match(.$Gene, gene_anotation$Probe_Id)])

    topde_tmp <-
      topde %>% dplyr::select(EntrezID, P.Value) %>% na.omit()

    topde_p <- tapply(topde_tmp$P.Value, topde_tmp$EntrezID, min)

    ranklist <- vector(mode = "numeric",
                       length = length(topde_p))
    names(ranklist) <- names(topde_p)

    # Temporary data frame with all values
    tmp <-
      filter(topde, EntrezID %in% names(topde_p) &
               P.Value %in% topde_p)

    # Adding values to the ranklist vector
    t_tmp <- tmp %>% dplyr::select(EntrezID, t)
    ranklist[as.character(t_tmp$EntrezID)] <- t_tmp$t

    # Initializing the genelist variable
    genelist <- tmp %>% dplyr::select(Gene, EntrezID, t, Name)
    rownames(genelist) <- as.character(1:nrow(genelist))

    # Sorting the ranklist
    ranklist <- sort(ranklist)

    # return(names(ranklist))

    # Geneset reactome
    # Muteing warning messages for now

    suppressMessages({
      suppressWarnings({
        # geneset_reactome
        geneset_reactome <- reactomePathways(names(ranklist))
        geneset_reactome <-
          geneset_reactome[intersect(names(geneset_reactome), toppw)]

        # Set seed
        set.seed(seed)

        # fgsea function
        fgseaRes <- fgsea(
          pathways = geneset_reactome,
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
        edge_entrez <-
          AnnotationDbi::select(
            org.Hs.eg.db,
            keys = edge,
            columns = "ENTREZID",
            keytype = "SYMBOL"
          )
        probe <-
          genelist$Gene[genelist$EntrezID %in% edge_entrez$ENTREZID]
        dat_selected <- transcriptomics_data[probe, , drop = FALSE]

        # Saving the data in a new frame
        saveRDS(dat_selected, "transcriptomics_selected.rds")

      })
    })

    return("Feature selection transcriptomics done!")

  }
