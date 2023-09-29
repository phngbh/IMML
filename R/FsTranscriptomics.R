#' Feature Selection Transcriptomics Data
#'
#' @description Feature selection for transcriptomics data using GSEA.
#'
#' @param trainIDs Set of training IDs from the `DataPartitioning()` function
#'   for the feature selection.
#' @param testIDs Set of testing IDs from the `DataPartitioning()` function for
#'   the feature selection.
#' @param dataIDs A data.frame with samples as rows and the data modalities as
#'   columns. It holds the data IDs of a sample for each modality. If for a
#'   sample there is no data for a modality, it has to be indicated by NA.
#' @param phenotypeIDs A data.frame with samples as rows and sample IDs as row
#'   names. Columns are phenotypes of interest. If for a sample no information
#'   about a phenotype is available, it has to be indicated by NA.
#' @param phenotype The name of the column in phenotypeIDs, which will be used
#'   in the analysis.
#' @param transcriptomicsData A data.frame with genes as rows and gene IDs as
#'   row names. Columns are samples with the sample IDs as column names.
#' @param geneAnnotation A data.frame of gene IDs and corresponding gene
#'   symbols, has at least 3 columns: "ID", "Symbol" and "EntrezID".
#' @param pathwayList A list of pathways used in GSEA.
#' @param seed The seed used for random number generation. Using the same seed
#'   ensures reproducibility.
#'
#' @return Returns transcriptomicsData subset to the selected features.
#'
#' @author Ulrich Asemann & Wilhelm Glaas
#'
#' @export

FsTranscriptomics <-
  function(trainIDs,
           testIDs,
           dataIDs,
           phenotypeIDs,
           phenotype,
           transcriptomicsData,
           geneAnnotation,
           pathwayList,
           seed = 123) {

    # Save row names as a column
    dataIDs <- dataIDs %>%
      tibble::rownames_to_column(var = "sampleIDs")
    phenotypeIDs <- phenotypeIDs %>%
      tibble::rownames_to_column(var = "sampleIDs")

    # Making a data frame with the phenotype info and transcriptomic IDs
    merged <- merge(dataIDs, phenotypeIDs, by = "sampleIDs")
    info <- merged %>% dplyr::select("Transcriptomics", any_of(phenotype)) %>%
      drop_na() %>% tibble::column_to_rownames(var = "Transcriptomics")

    # Removing all samples not relevant for transcriptomics feature selection
    modelIDs <- dataIDs %>% drop_na() %>% dplyr::select("Transcriptomics")
    info <- info %>% filter(!(row.names(info) %in% modelIDs$Transcriptomics))

    # Creating a model matrix
    target <- info %>% transmute(inc3 = as.character(inc3))

    # some of transcIDs which are in dataIDs are not in transcriptomicsData
    data <- transcriptomicsData[, rownames(target)[rownames(target) %in% colnames(transcriptomicsData)]]
    rm(transcriptomicsData)

    modmatrix <- model.matrix( ~ 0 + ., data = target)
    modmatrix <- modmatrix[colnames(data), ]

    gsea <- list(edge = list(), ilmn = list(), auc = list(), pathway = list())
    set.seed(seed)


    for (i in 1:length(trainIDs)){
      cat("Iter ",i,"\n")

      # some of transcIDs which are in trainIDs are not in data
      trainIDs[[i]] <- trainIDs[[i]][!(trainIDs[[i]] %in% setdiff(trainIDs[[i]], colnames(data)))]
      testIDs[[i]] <- testIDs[[i]][!(testIDs[[i]] %in% setdiff(testIDs[[i]], colnames(data)))]

      dat_tmp <- data[, as.character(trainIDs[[i]])]

      cat("...DE analysis\n")
      mod_tmp <- modmatrix[as.character(trainIDs[[i]]),]
      if( any(table(mod_tmp[,1]) < 2) | any(table(mod_tmp[,2]) < 2)){
        cat("......One of the factor levels has less than 2 observations => Stop!")
        next
      }
      fit <- lmFit(dat_tmp, mod_tmp)
      contrast <- makeContrasts(contrasts = paste0(phenotype,"1-",phenotype,"0"), levels = colnames(coef(fit)))
      tmp <- contrasts.fit(fit, contrast)
      tmp <- eBayes(tmp)
      topde <- topTable(tmp, sort.by = "P", n = Inf) %>% mutate(ID = rownames(.)) %>%
        mutate(Name = geneAnnotation$Symbol[match(.$ID,geneAnnotation$ID)],
               EntrezID = geneAnnotation$EntrezID[match(.$ID,geneAnnotation$ID)])

      cat("...GSEA\n")
      topde.tmp <- dplyr::select(topde, EntrezID, P.Value) %>% na.omit()
      topde.p <- tapply(topde.tmp$P.Value, topde.tmp$EntrezID, min) # Only select the most significant probe for each gene
      ranklist <- vector(mode = "numeric", length = length(topde.p))
      names(ranklist) <- names(topde.p)
      ilmnlist <- vector(mode = "character", length = length(topde.p))
      for(l in 1:length(topde.p)){
        t <- filter(topde, EntrezID == names(topde.p)[l] & P.Value == topde.p[l])$t
        il <- filter(topde, EntrezID == names(topde.p)[l] & P.Value == topde.p[l])$ID
        ranklist[l] = t
        ilmnlist[l] = il
      }
      genelist_tmp <- data.frame(Entrez = names(ranklist), Probe = ilmnlist, t = ranklist) %>%
        mutate(Name = topde$Name[match(.$Entrez, topde$EntrezID)])
      ranklist <- sort(ranklist)
      geneset_reactome <- reactomePathways(names(ranklist))
      geneset_reactome <- geneset_reactome[intersect(names(geneset_reactome),pathwayList)]
      set.seed(seed)
      fgseaRes_tmp <- fgsea(pathways = geneset_reactome,
                            stats    = ranklist,
                            minSize  = 15,
                            maxSize  = 200,
                            eps = 0) %>% arrange(pval) %>% filter(padj < 0.1)
      if (nrow(fgseaRes_tmp) == 0){
        message("No significant pathway found")
        next
      }
      edge_tmp = fgseaRes_tmp$leadingEdge %>% unlist() %>% unique()
      gsea$edge[[i]] = edge_tmp
      gsea$ilmn[[i]] = genelist_tmp$Probe[match(edge_tmp, genelist_tmp$Entrez)]
      gsea$pathway[[i]] = fgseaRes_tmp

      cat("...Elastic net\n")

      if(is.null(gsea$ilmn[[i]])){
        message("No significant pathway found")
        next
      }

      probelistTmp = gsea$ilmn[[i]]
      x_train = data[probelistTmp, as.character(trainIDs[[i]])] %>% t()
      y_train = target[as.character(trainIDs[[i]]), phenotype]
      y_train = ifelse(y_train == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
      x_test = data[probelistTmp, as.character(testIDs[[i]])] %>% t()
      y_test = target[as.character(testIDs[[i]]), phenotype]
      y_test = ifelse(y_test == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))

      my_control <- trainControl(
        method="repeatedcv",
        number=5,
        repeats = 2,
        savePredictions="final",
        classProbs=TRUE,
        summaryFunction=twoClassSummary,
        sampling = "smote",
        allowParallel = T
      )
      set.seed(seed)
      fit <- caret::train(x = x_train,
                          y = y_train,
                          method="glmnet",
                          metric="ROC",
                          tuneLength = 20,
                          maximize = T,
                          trControl=my_control,
                          importance = TRUE)
      pred = predict(fit, x_test, s = "lambda.min", type = "prob")$One
      roc <- roc(response = y_test, predictor = pred, levels = c("Zero","One"))
      auc = auc(roc)
      gsea$auc[[i]] = auc

    }

    #Select the top significant pathways
    keep = vector("logical",length = length(gsea$auc))
    for (i in 1:length(keep)){
      if (is.null(gsea$auc[[i]])){
        next
      } else if (gsea$auc[[i]] <= 0.5){
        next
      } else {
        keep[i] = TRUE
      }
    }

    pwlist = gsea$pathway[keep]
    pathways = list(pw = list(), pval = list(), NES = list())
    for (i in 1: length(pwlist)){
      if(is.null(pwlist[[i]])){
        next
      }
      pathways$pw[[i]] = pwlist[[i]]$pathway
      pathways$pval[[i]] = pwlist[[i]]$pval
      pathways$NES[[i]] = pwlist[[i]]$NES
    }
    pwlist_agg = aggregateRanks(pathways$pw)
    pwlist_agg$adjP = pwlist_agg$Score*length(pathways$pw)
    pwlist_agg$adjP = p.adjust(pwlist_agg$adjP, method = "fdr")
    toppw = rownames(filter(pwlist_agg, adjP < 0.05))


    #Final gene set enrichment analysis on the selected pathways
    cat("Final GSEA\n")
    fit = lmFit(data, modmatrix)
    contrast = makeContrasts(contrasts = paste0(phenotype,"1-",phenotype,"0"), levels = colnames(coef(fit)))
    tmp <- contrasts.fit(fit, contrast)
    tmp <- eBayes(tmp)
    topde <- topTable(tmp, sort.by = "P", n = Inf)
    topde = topde %>% mutate(ID = rownames(topde)) %>%
      mutate(Name = geneAnnotation$Symbol[match(.$ID,geneAnnotation$ID)],
             EntrezID = geneAnnotation$EntrezID[match(.$ID,geneAnnotation$ID)])
    topde.tmp = dplyr::select(topde, EntrezID, P.Value) %>% na.omit()
    tmp = tapply(topde.tmp$P.Value, topde.tmp$EntrezID, min)
    ranklist = vector(mode = "numeric", length = length(tmp))
    names(ranklist) = names(tmp)
    ilmnlist = vector(mode = "character", length = length(tmp))
    for(i in 1:length(tmp)){
      t = filter(topde, EntrezID == names(tmp)[i] & P.Value == tmp[i])$t
      il = filter(topde, EntrezID == names(tmp)[i] & P.Value == tmp[i])$ID
      ranklist[i] = t
      ilmnlist[i] = il
    }
    genelist = data.frame(Entrez = names(ranklist), Probe = ilmnlist, t = ranklist)
    ranklist = sort(ranklist)
    geneset_reactome = reactomePathways(names(ranklist))
    geneset_reactome = geneset_reactome[intersect(names(geneset_reactome),toppw)]

    set.seed(seed)
    fgseaRes <- fgsea(pathways = geneset_reactome,
                      stats    = ranklist,
                      minSize  = 15,
                      maxSize  = 200) %>% arrange(pval)

    #Extract selected features for training
    edge = fgseaRes$leadingEdge %>% unlist() %>% unique()
    #edge_entrez = AnnotationDbi::select(org.Hs.eg.db, keys=edge, columns="ENTREZID", keytype="SYMBOL")
    probe = genelist$Probe[genelist$Entrez %in% edge]

    return(list(selected_feature = probe, gsea_result = fgseaRes))


    # Saving the data in a new frame
    #saveRDS(dataSelected, "transcriptomics_selected.rds")

    return("Feature selection transcriptomics done!")
  }
