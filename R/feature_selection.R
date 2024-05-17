suppressMessages(library(tibble))
#suppressMessages(library(tidyverse))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(xtable))
suppressMessages(library(stringr))
suppressMessages(library(scales))
suppressMessages(library(devtools))
#suppressMessages(library(reactome.db))
suppressMessages(library(fgsea))
suppressMessages(library(pROC))
suppressMessages(library(RobustRankAggreg))
suppressMessages(library(globaltest))
suppressMessages(library(RaMP))
suppressMessages(library(caret))
suppressMessages(library(limma))

data_partitioning <-
  function(phenotypeIDs,
           dataIDs,
           phenotype,
           partitioning = 0.8,
           k = 5,
           iter = 100,
           seed = 123) {
    # Setting seed
    set.seed(seed)
    
    # Defining variables
    trainingIDs <- list()
    trainSingleFeatureList <- list()
    testSingleFeatureList <- list()
    testFeatureList <- list()
    trainFeatureList <- list()
    finalList <- list()
    
    # List of feature selection modalities
    listNames <- colnames(dataIDs)
    
    # Selecting IDs used for training and feature selection
    IDs <-
      phenotypeIDs %>%
      dplyr::select(any_of(phenotype)) %>%
      drop_na() %>%
      rownames() %>%
      as.double()
    
    # Save row names as a column
    dataIDs <- dataIDs %>% tibble::rownames_to_column(var = "sampleID")
    
    # Data frame with the selected IDs
    newSampleIDs <- dataIDs %>% filter(sampleID %in% IDs)
    
    # Building the model training IDs list
    
    # Selecting sample IDs with IDs for each modality
    trainingSampleIDs <-
      newSampleIDs %>% drop_na() %>% dplyr::select(sampleID)
    
    # Building a list with the modeltraining IDs
    trainingIDs <- as.double(trainingSampleIDs$sampleID)
    
    iterations <-
      caret::createDataPartition(trainingIDs, times = iter, p = partitioning)
    
    returnList <- list()
    for (i in 1:iter) {
      tmpFoldList <- list()
      
      testIDs <- trainingIDs[-iterations[[i]]]
      trainIDs <- createFolds(trainingIDs[iterations[[i]]], k = k)
      trainIDs <- lapply(trainIDs, function(x) trainingIDs[x])
      
      tmpFoldList[['Training']] <- trainIDs
      tmpFoldList[['Testing']] <- testIDs
      
      returnList[[paste0('Iteration', i)]] <- tmpFoldList
    }
    # append to the final list
    finalList <- c(finalList, "modelTraining" = list(returnList))
    
    # Remove the modeltraining IDs from the newSampleIDs
    newSampleIDs <-
      newSampleIDs %>% filter(!(sampleID %in% trainingIDs))
    
    fsList = list()
    # Building a feature selection list for each modality
    for (i in 1:length(listNames)) {
      # Select one feature and remove all NA values from the newSampleIDs list
      frame <-
        newSampleIDs %>% dplyr::select(listNames[i]) %>% drop_na()
      
      fsVec = frame[[1]]
      fsList[[listNames[i]]] = fsVec
    }
    finalList[['featureSelection']] = fsList
    # Return the completed list
    return(finalList)
  }

standardise = function(
    ###Scale datasets into standard distribution
  df = NULL #dataframe or matrix to be transformed
){
  rm = apply(as.data.frame(df), 2, function(x) sd(x) == 0)
  df = df[,!rm]
  df = apply(as.data.frame(df), 2, function(x) (x - mean(x))/sd(x)) %>% as.matrix()
  return(df)
}

get_pathways <- function(
    # Get all pathways in RaMP database that are associated with the analytes of interest
  # Output: A list of pathways and their member analytes (with local dataset IDs)
  database_id = NULL, # character vector of database IDs of the analytes, in form '<database>:<id>', eg: kegg:xxxxx, entrez:xxxxx, uniprot:xxxx...
  local_id = NULL # character vector of the same length as database_id providing corresponding local data ID of the respective analytes
){
  # Quick check
  if(length(database_id) != length(local_id)){
    stop("Number of database IDs is different from number of local IDs. Please check again")
  }
  
  # Set up your connection to the RaMP2.0 database:
  ramp <- RaMP()
  pathwaydf <- getPathwayFromAnalyte(ramp, analytes = database_id) %>% unite(col = "pathwayName", pathwaySource, pathwayName, sep = ":")
  pathwaylist <- list()
  for (i in unique(pathwaydf$pathwayName)){
    .df <- filter(pathwaydf, pathwayName == i)
    if (nrow(.df) < 5 | length(unique(.df$inputId)) < 5){
      next
    }
    lid <- local_id[match(.df$inputId, database_id)]
    pathwaylist[[i]] <- unique(lid)
  }
  
  return(pathwaylist)
}

fit_elnet <- function(
    # Fit an elastic net model and output prediction performance on a test set
  # Output: test AUROC value
  data = NULL, # data matrix of shape (n_samples, n_features)
  target = NULL, # vector of target variable (0 = ctrl, 1 = case)
  train_samples = NULL, # character vector of training samples IDs
  variables = NULL, # character vector of variables to use
  seed = NULL
){
  
  x_train = data[train_samples,variables, drop = F]
  if (ncol(x_train) == 1){
    x_train <- cbind(x_train,ranv = 1)
  }
  train_ind <- match(train_samples, rownames(data))
  y_train = target[train_ind]
  y_train = ifelse(y_train == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
  test_samples = setdiff(rownames(data), train_samples)
  x_test = data[test_samples,variables, drop = F] 
  if (ncol(x_test) == 1){
    x_test <- cbind(x_test, ranv = 1)
  }
  test_ind <- setdiff(1:nrow(data), train_ind)
  y_test = target[test_ind]
  y_test = ifelse(y_test == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
  
  LogLoss <- function(pred, true, eps = 1e-15, weights = NULL) {
    pred = pmin(pmax(pred, eps), 1 - eps) # Bound the results
    if (is.null(weights)) {
      return(-(sum(
        true * log(pred) + (1 - true) * log(1 - pred)
      )) / length(true))
    } else{
      return(-weighted.mean(true * log(pred) + (1 - true) * log(1 - pred), weights))
    }
  }
  
  caretLogLoss <- function(data, lev = NULL, model = NULL) {
    cls <- levels(data$obs) #find class names
    loss <- LogLoss(
      pred = data[, cls[2]],
      true = as.numeric(data$obs) - 1,
      weights = data$weights
    )
    names(loss) <- c('myLogLoss')
    loss
  }
  
  sampling = NULL
  sumFunc = caretLogLoss
  metric = "myLogLoss"
  maximize = F
  
  my_control <- trainControl(
    method="repeatedcv",
    number=5,
    repeats = 4,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=sumFunc,
    sampling = sampling,
    allowParallel = T
  )
  weights = ifelse(y_train == "One", table(y_train)[[2]]/table(y_train)[[1]], 1)
  
  if (any(rownames(x_train) != names(y_train)) | any(rownames(x_test) != names(y_test))){
    message("Samples in train and test sets do not match!\n")
    next
  }
  set.seed(seed)
  fit <- caret::train(x = x_train,
                      y = y_train,
                      method="glmnet", 
                      metric=metric,
                      tuneLength = 20,
                      weights = weights,
                      maximize = maximize,
                      trControl=my_control,
                      importance = TRUE)
  pred = predict(fit, x_test, s = "lambda.min", type = "prob")$One
  roc <- roc(response = y_test, predictor = pred, levels = c("Zero","One"))
  auc = auc(roc)
  
  return(auc)
}

featureSelection_untargeted <- function(
    # Do feature selection for untargeted transcriptomics, proteomics and metabolomics data
  # Output: a vector of selected features
  data = NULL, # matrix of processed data
  target = NULL, # dataframe of target variable and (optional) technical variables
  target_name = NULL, # string of target name (must be a column name in target dataframe)
  #feature_annotation = NULL, # dataframe of feature IDs and corresponding gene symbols, has at least 3 columns: "ID", "Symbol" and "EntrezID"/"CheBI"
  pathway_list = NULL, # character vector of pathways to test
  seed = 993, # random seed
  resampling = TRUE, # Whether should do the analysis in many resamples
  n_iterations = 100, # number of resamples (if resampling == TRUE)
  p = 0.8,
  GSEA_FDR = 0.2
){
  
  # Make sure the dimensions are correct
  if (any(colnames(data) != rownames(target))){
    cat("Columns of data do not match rows of target dataframe => Adjusting the tables\n")
    .samples <- intersect(colnames(data), rownames(target))
    data <- data[, .samples]
    target <- target[.samples,, drop = F]
  } 
  
  #Do gene set enrichment analysis
  #modmatrix <- model.matrix(~ 0 + ., data=target)
  target[,target_name] <- as.factor(target[,target_name])
  f <- reformulate(termlabels = colnames(target), intercept = F)
  modmatrix <- model.matrix(f, data = target)
  
  if (resampling){
    cat("Resampling approach is set => Do feature selection on resamples\n")
    
    if (ncol(data) < 100){
      message("Number of samples is less than 100 so please consider setting lower number of resamples or do not use resampling at all.\n")
    }
    # Initiate result list
    res <- list(auc_gsea = list(), auc_de = list(), auc_union = list(), pathway = list(), gsea_probe = list(), gsea_probe_rra = list(), d_probe = list(), d_probe_rra = list(), resampling = T)
    
    # Make resamples for feature selection
    set.seed(seed)
    resamples <- createDataPartition(y = target[,target_name], times = n_iterations, p = p)
    for (i in 1:length(resamples)){
      cat("Iter ",i,"\n")
      data_tmp <- data[,resamples[[i]]]
      
      cat("...DE analysis\n")
      mod_tmp <- modmatrix[resamples[[i]],,drop = F]
      if( any(table(mod_tmp[,1]) < 2) | any(table(mod_tmp[,2]) < 2)){
        cat("......One of the factor levels has less than 2 observations => Stop!\n")
        next
      }
      fit <- lmFit(data_tmp, mod_tmp)
      ctrst <- paste0(target_name,"1 - ",target_name,"0")
      contrast <- makeContrasts(contrasts = ctrst, levels = colnames(coef(fit)))
      tmp <- contrasts.fit(fit, contrast)
      tmp <- eBayes(tmp)
      topde <- topTable(tmp, sort.by = "P", n = Inf) %>% mutate(ID = rownames(.))
      res$d_probe[[i]] <- topde$ID[topde$adj.P.Val < 0.05]
      cat(paste0("......There are ",length(res$d_probe[[i]])," significant probes\n"))
      
      cat("...GSEA\n")
      ranklist <- topde$t
      names(ranklist) <- topde$ID
      ranklist <- sort(ranklist)
      set.seed(seed)
      fgseaRes_tmp <- fgsea(pathways = pathway_list,
                            stats    = ranklist,
                            minSize  = 5,
                            maxSize  = 200,
                            eps = 0) %>% arrange(pval) %>% filter(padj < GSEA_FDR)
      cat(paste0("......There are ",nrow(fgseaRes_tmp)," significant pathways\n"))
      res$pathway[[i]] <- fgseaRes_tmp$pathway
      edge_tmp <- fgseaRes_tmp$leadingEdge %>% unlist() %>% unique()
      cat(paste0("......There are ",length(edge_tmp)," leading edge probes associated with the pathways\n"))
      res$gsea_probe[[i]] <- edge_tmp
      
      cat("...Elastic net\n")
      if (length(res$d_probe[[i]]) > 0){
        cat("......Evaluate prediction performance for differentially expressed probes\n")
        auc <- fit_elnet(data = t(data), train_samples = colnames(data)[resamples[[i]]], target = target[,target_name], variables = res$d_probe[[i]], seed = seed)
        res$auc_de[[i]] = auc
        cat(paste0("......Test AUROC is ",auc,"\n"))
      } else {
        cat("......There is no differentially expressed probes to analyse\n")
      }
      if (length(edge_tmp) > 0){
        cat("......Evaluate prediction performance for GSEA-derived probes\n")
        auc <- fit_elnet(data = t(data), train_samples = colnames(data)[resamples[[i]]], target = target[,target_name], variables = edge_tmp, seed = seed)
        res$auc_gsea[[i]] = auc
        cat(paste0("......Test AUROC is ",auc,"\n"))
      } else {
        cat("......There is no leading edge probes to analyse\n")
      }
      # Evaluate the combined set of selected features
      if (length(edge_tmp) > 0 & length(res$d_probe[[i]]) > 0){
        cat("......Evaluate prediction performance of the union of both types of features\n")
        auc <- fit_elnet(data = t(data), train_samples = colnames(data)[resamples[[i]]], target = target[,target_name], variables = union(edge_tmp, res$d_probe[[i]]), seed = seed)
        res$auc_union[[i]] = auc
        cat(paste0("......Test AUROC is ",auc,"\n"))
      }
      
    }
    
    #Aggregate significant probes
    cat("Aggregate lists of differentially expressed probes to extract final significant probes\n")
    if (length(res$auc_de) == 0){
      cat("...There is no significant probe list to aggregate\n")
      res$d_probe_rra <- NA
    } else {
      #Select the iterations with good prediction performance
      keep = vector("logical",length = length(res$auc_de))
      for (i in 1:length(keep)){
        if (is.null(res$auc_de[[i]])){
          next
        } else if (res$auc_de[[i]] <= 0.5){
          next
        } else {
          keep[i] = TRUE
        }
      }
      probelist <- res$d_probe[keep]
      if (length(probelist) == 0){
        cat("...There is no significant probe list to aggregate\n")
        res$d_probe_rra <- NA
      } else {
        set.seed(seed)
        de_probes_agg = aggregateRanks(probelist)
        de_probes_agg$adjP = de_probes_agg$Score*length(probelist)
        de_probes_agg$adjP = p.adjust(de_probes_agg$adjP, method = "fdr")
        aggregated_probes = rownames(filter(de_probes_agg, adjP < 0.05))
        cat(paste0("...There are ",length(aggregated_probes)," aggregated probes\n"))
        res$d_probe_rra <- aggregated_probes
      }
    }
    
    cat("Aggregate lists of GSEA-derived pathways to extract final significant probes\n")
    if (length(res$auc_gsea) == 0){
      cat("...There is no significant pathway list to aggregate\n")
      res$gsea_probe_rra <- NA
    } else {
      #Select the top significant pathways
      keep = vector("logical",length = length(res$auc_gsea))
      for (i in 1:length(keep)){
        if (is.null(res$auc_gsea[[i]])){
          next
        } else if (res$auc_gsea[[i]] <= 0.5){
          next
        } else {
          keep[i] = TRUE
        }
      }
      pwlist = res$pathway[keep]
      if (length(pwlist) == 0){
        cat("...There is no significant pathway list to aggregate\n")
        res$gsea_probe_rra <- NA
      } else {
        set.seed(seed)
        pwlist_agg = aggregateRanks(pwlist)
        pwlist_agg$adjP = pwlist_agg$Score*length(pwlist)
        pwlist_agg$adjP = p.adjust(pwlist_agg$adjP, method = "fdr")
        toppw = rownames(filter(pwlist_agg, adjP < 0.05))
        cat(paste0("...There are ",length(toppw)," aggregated pathways\n"))
        if (length(toppw) == 0){
          res$gsea_probe_rra <- NA
        } else {
          #Final gene set enrichment analysis on the selected pathways
          cat("...Final GSEA to extract final leading edge probes\n")
          fit = lmFit(data, modmatrix)
          contrast = makeContrasts(contrasts = paste0(target_name,"1-",target_name,"0"), levels = colnames(coef(fit)))
          tmp <- contrasts.fit(fit, contrast)
          tmp <- eBayes(tmp)
          topde <- topTable(tmp, sort.by = "P", n = Inf) %>% mutate(ID = rownames(.))
          ranklist <- topde$t
          names(ranklist) <- topde$ID
          ranklist <- sort(ranklist)
          set.seed(seed)
          fgseaRes <- fgsea(pathways = pathway_list[toppw],
                            stats    = ranklist,
                            minSize  = 5,
                            maxSize  = 200) %>% arrange(pval) 
          edge = fgseaRes$leadingEdge %>% unlist() %>% unique()
          cat(paste0("......There are ",length(edge)," final leading edge probes\n"))
          res$gsea_probe_rra <- edge
        }
      }
      
    }
    return(res)
    
  } else {
    cat("Resampling approach is not set => Do feature selection on whole feature selection set\n")
    res <- list(pathway = list(), gsea_probe = list(),  d_probe = list(), resampling = F)
    
    cat("DE analysis\n")
    fit = lmFit(data, modmatrix)
    contrast = makeContrasts(contrasts = paste0(target_name,"1-",target_name,"0"), levels = colnames(coef(fit)))
    tmp <- contrasts.fit(fit, contrast)
    tmp <- eBayes(tmp)
    topde <- topTable(tmp, sort.by = "P", n = Inf) %>% mutate(ID = rownames(.))
    res$d_probe <- topde$ID[topde$adj.P.Val < 0.05]
    cat(paste0("......There are ",length(res$d_probe)," significant probes\n"))
    
    cat("GSEA\n")
    ranklist <- topde$t
    names(ranklist) <- topde$ID
    ranklist <- sort(ranklist)
    set.seed(seed)
    fgseaRes <- fgsea(pathways = pathway_list,
                      stats    = ranklist,
                      minSize  = 5,
                      maxSize  = 200) %>% arrange(pval) %>% filter(padj < 0.05)
    cat(paste0("......There are ",nrow(fgseaRes)," significant pathways\n"))
    res$pathway <- fgseaRes$pathway
    edge = fgseaRes$leadingEdge %>% unlist() %>% unique() 
    cat(paste0("......There are ",length(edge)," final leading edge probes\n"))
    res$gsea_probe <- edge
    
    return(res)
  }
  
}

featureSelection_targeted <- function(
    # Do feature selection for targeted transcriptomics, proteomics and metabolomics data
  # Output: a vector of selected features
  data = NULL, # matrix of processed data
  target = NULL, # dataframe of target variable and (optional) technical variables
  target_name = NULL, # string of target name (must be a column name in target dataframe)
  pathway_list = NULL, # list of signaling pathway and their analyte members
  seed = 993, # random seed
  resampling = TRUE, # Whether should do the analysis in many resamples
  n_iterations = 100, # number of resamples (if resampling == TRUE)
  p = 0.8,
  MSEA_FDR = 0.2 # FDR cut off for MSEA in each resample
){
  
  # Make sure the dimensions are correct
  if (any(colnames(data) != rownames(target))){
    cat("Columns of data do not match rows of target dataframe => Adjusting the tables\n")
    .samples <- intersect(colnames(data), rownames(target))
    data <- data[, .samples]
    target <- target[.samples,, drop = F]
  } 
  
  # Make model matrix for univariate analysis
  target[,target_name] <- as.factor(target[,target_name])
  f <- reformulate(termlabels = colnames(target), intercept = F)
  modmatrix <- model.matrix(f, data = target)
  
  if (resampling){
    cat("Resampling approach is set => Do feature selection on resamples\n")
    
    if (ncol(data) < 100){
      message("Number of samples is less than 100 so please consider setting lower number of resamples or do not use resampling at all.\n")
    }
    
    # Initiate result list
    res <- list(auc_msea = list(), auc_de = list(), auc_union = list(), pathway = list(), msea_m = list(), msea_m_rra = list(), d_m = list(), d_m_rra = list(), resampling = T)
    
    # Make resamples of feature selection set
    set.seed(seed)
    resamples <- createDataPartition(y = target[,target_name], times = n_iterations, p = p)
    
    for (i in 1:length(resamples)){
      cat("Iter ",i,"\n")
      
      cat("...Differential expression analysis to get significant univariate probes\n")
      data_tmp <- data[,resamples[[i]]]
      mod_tmp <- modmatrix[resamples[[i]],,drop = F]
      if( any(table(mod_tmp[,1]) < 2) | any(table(mod_tmp[,2]) < 2)){
        cat("......One of the factor levels has less than 2 observations => Stop!\n")
        next
      }
      fit <- lmFit(data_tmp, mod_tmp)
      ctrst <- paste0(target_name,"1 - ",target_name,"0")
      contrast <- makeContrasts(contrasts = ctrst, levels = colnames(coef(fit)))
      tmp <- contrasts.fit(fit, contrast)
      tmp <- eBayes(tmp)
      topde <- topTable(tmp, sort.by = "P", n = Inf) %>% mutate(ID = rownames(.)) 
      
      res$d_m[[i]] <- topde$ID[topde$adj.P.Val < 0.05]
      cat(paste0("......There are ",length(res$d_m[[i]])," significant probes\n"))
      
      cat("...Pathway enrichment analysis using global test\n")
      data_tmp <- data[,resamples[[i]]] %>% t() 
      y <- target[resamples[[i]],target_name]
      cat("......Iterate through all pathways\n")
      pplist <- list()
      for (p in names(pathway_list)){
        mat <- data_tmp[,pathway_list[[p]]]
        res_gt <- gt(y, mat, model = "logistic")
        pval <- p.value(res_gt)
        pplist[[p]] <- pval
      }
      adjplist <- unlist(pplist) %>% p.adjust(method = "fdr")
      sig_mset <- adjplist[adjplist < MSEA_FDR] %>% sort() %>% names()
      cat(paste0("......There are ",length(sig_mset)," significant pathways\n"))
      res$pathway[[i]] = sig_mset
      sig_m <- pathway_list[sig_mset] %>% unlist() %>% unique() 
      cat(paste0("......There are ",length(sig_m)," unique probes within the pathways\n"))
      res$msea_m[[i]] <- sig_m
      
      cat("...Elastic net\n")
      if (length(res$d_m[[i]]) > 0){
        cat("......Evaluate prediction performance for differentially expressed probes\n")
        auc <- fit_elnet(data = t(data), train_samples = colnames(data)[resamples[[i]]], target = target[,target_name], variables = res$d_m[[i]], seed = seed)
        res$auc_de[[i]] = auc
        cat(paste0("......Test AUROC is ",auc,"\n"))
      } else {
        cat("......There is no differentially expressed probes to analyse\n")
      }
      if (length(sig_m) > 0){
        cat("......Evaluate prediction performance for global test-derived probes\n")
        auc <- fit_elnet(data = t(data), train_samples = colnames(data)[resamples[[i]]], target = target[,target_name], variables = sig_m, seed = seed)
        res$auc_msea[[i]] = auc
        cat(paste0("......Test AUROC is ",auc,"\n"))
      } else {
        cat("......There is no leading edge probes to analyse\n")
      }
      # Evaluate the combined set of selected features
      if (length(sig_m) > 0 & length(res$d_m[[i]]) > 0){
        cat("......Evaluate prediction performance of the union of both types of features\n")
        auc <- fit_elnet(data = t(data), train_samples = colnames(data)[resamples[[i]]], target = target[,target_name], variables = union(sig_m, res$d_m[[i]]), seed = seed)
        res$auc_union[[i]] = auc
        cat(paste0("......Test AUROC is ",auc,"\n"))
      }
      
    }
    
    #Aggregate significant probes
    cat("Aggregate lists of differentially expressed probes to extract final significant probes\n")
    if (length(res$auc_de) == 0){
      cat("...There is no significant probe list to aggregate\n")
      res$d_m_rra <- NA
    } else {
      #Select the iterations with good prediction performance
      keep = vector("logical",length = length(res$auc_de))
      for (i in 1:length(keep)){
        if (is.null(res$auc_de[[i]])){
          next
        } else if (res$auc_de[[i]] <= 0.5){
          next
        } else {
          keep[i] = TRUE
        }
      }
      probelist <- res$d_m[keep]
      if (length(probelist) == 0){
        cat("...There is no significant probe list to aggregate\n")
        res$d_m_rra <- NA
      } else {
        set.seed(seed)
        de_probes_agg = aggregateRanks(probelist)
        de_probes_agg$adjP = de_probes_agg$Score*length(probelist)
        de_probes_agg$adjP = p.adjust(de_probes_agg$adjP, method = "fdr")
        aggregated_probes = rownames(filter(de_probes_agg, adjP < 0.05))
        cat(paste0("...There are ",length(aggregated_probes)," aggregated probes\n"))
        res$d_m_rra <- aggregated_probes
      }
    }
    cat("Aggregate lists of global test-derived pathways to extract final significant probes\n")
    if (length(res$auc_msea) == 0){
      cat("...There is no significant pathway list to aggregate\n")
      res$msea_m_rra <- NA
    } else {
      #Select the top significant pathways
      keep = vector("logical",length = length(res$auc_msea))
      for (i in 1:length(keep)){
        if (is.null(res$auc_msea[[i]])){
          next
        } else if (res$auc_msea[[i]] <= 0.5){
          next
        } else {
          keep[i] = TRUE
        }
      }
      pwlist = res$pathway[keep]
      if (length(pwlist) == 0){
        cat("...There is no significant pathway list to aggregate\n")
        res$msea_m_rra <- NA
      } else {
        set.seed(seed)
        pwlist_agg = aggregateRanks(pwlist)
        pwlist_agg$adjP = pwlist_agg$Score*length(pwlist)
        pwlist_agg$adjP = p.adjust(pwlist_agg$adjP, method = "fdr")
        toppw = rownames(filter(pwlist_agg, adjP < 0.05))
        cat(paste0("...There are ",length(toppw)," aggregated pathways\n"))
        if (length(toppw) == 0){
          res$msea_m_rra <- NA
        } else {
          sig_m <- pathway_list[toppw] %>% unlist() %>% unique() 
          cat(paste0("...There are ",length(sig_m)," unique probes within the pathways\n"))
          res$msea_m_rra <- sig_m
        }
      }
      
    }
    
    return(res)
    
  } else {
    cat("Resampling approach is not set => Do feature selection on the whole dataset\n")
    
    # Initiate result list
    res <- list()
    
    cat("Differential expression analysis to get significant univariate probes\n")
    data_tmp <- data
    mod_tmp <- modmatrix
    if( any(table(mod_tmp[,1]) < 2) | any(table(mod_tmp[,2]) < 2)){
      cat("...One of the factor levels has less than 2 observations => Stop!\n")
      next
    }
    fit <- lmFit(data_tmp, mod_tmp)
    ctrst <- paste0(target_name,"1 - ",target_name,"0")
    contrast <- makeContrasts(contrasts = ctrst, levels = colnames(coef(fit)))
    tmp <- contrasts.fit(fit, contrast)
    tmp <- eBayes(tmp)
    topde <- topTable(tmp, sort.by = "P", n = Inf) %>% mutate(ID = rownames(.)) 
    
    res[["d_m"]] <- topde$ID[topde$adj.P.Val < 0.05]
    cat(paste0("...There are ",length(res[["d_m"]])," significant probes\n"))
    
    cat("Pathway enrichment analysis using global test\n")
    data_tmp <- data %>% t() 
    y <- target[,target_name]
    cat("...Iterate through all pathways\n")
    pplist <- list()
    for (p in names(pathway_list)){
      mat <- data_tmp[,pathway_list[[p]]]
      res_gt <- gt(y, mat, model = "logistic")
      pval <- p.value(res_gt)
      pplist[[p]] <- pval
    }
    adjplist <- unlist(pplist) %>% p.adjust(method = "fdr")
    sig_mset <- adjplist[adjplist < MSEA_FDR] %>% sort() %>% names()
    cat(paste0("...There are ",length(sig_mset)," significant pathways\n"))
    res[["pathway"]] = sig_mset
    sig_m <- pathway_list[sig_mset] %>% unlist() %>% unique() 
    cat(paste0("......There are ",length(sig_m)," unique probes within the pathways\n"))
    res[["msea_m"]] <- sig_m
    res[["resampling"]] = F
    
    return(res)
  }
  
}

featureSelection_clinical <- function(
    # Do feature selection for clinical data
  # Output: a vector of selected features
  data = NULL, # matrix of processed data
  target = NULL, # dataframe of target variable and (optional) technical variables
  target_name = NULL, # string of target name (must be a column name in target dataframe)
  to_remove = NULL, # character vector of variables to be removed from analysis (usually ones that were used to compute the target)
  seed = 993, # random seed
  resampling = TRUE, # Whether should do the analysis in many resamples
  n_iterations = 100, # number of resamples (if resampling == TRUE)
  p = 0.8
){
  
  # Make sure the dimensions are correct
  if (any(rownames(data) != rownames(target))){
    cat("Columns of data do not match rows of target dataframe => Adjusting the tables\n")
    .samples <- intersect(rownames(data),rownames(target))
    data <- data[.samples,]
    target <- target[.samples,,drop = F]
  } 
  
  # Remove target-related variables
  if (!is.null(to_remove)){
    cat("Remove target-related variables\n")
    data <- dplyr::select(data, -all_of(to_remove))
  }
  
  cat("Prepare resamples\n")
  target[,target_name] <- ifelse(target[,target_name] == 1, "One", "Zero") %>% factor(levels = c("One","Zero"))
  resamples <- createDataPartition(y = target[,target_name], times = n_iterations, p = p)
  
  cat("Set up training parameters\n")
  LogLoss <- function(pred, true, eps = 1e-15, weights = NULL) {
    pred = pmin(pmax(pred, eps), 1 - eps) # Bound the results
    if (is.null(weights)) {
      return(-(sum(
        true * log(pred) + (1 - true) * log(1 - pred)
      )) / length(true))
    } else{
      return(-weighted.mean(true * log(pred) + (1 - true) * log(1 - pred), weights))
    }
  }
  
  caretLogLoss <- function(data, lev = NULL, model = NULL) {
    cls <- levels(data$obs) #find class names
    loss <- LogLoss(
      pred = data[, cls[2]],
      true = as.numeric(data$obs) - 1,
      weights = data$weights
    )
    names(loss) <- c('myLogLoss')
    loss
  }
  
  sampling = NULL
  sumFunc = caretLogLoss
  metric = "myLogLoss"
  maximize = F
  
  my_control <- trainControl(
    method="repeatedcv",
    number=5,
    repeats = 4,
    savePredictions="final",
    classProbs=TRUE,
    summaryFunction=sumFunc,
    sampling = sampling,
    allowParallel = T
  )
  
  perf_list = list()
  var_list = list()
  for (i in 1:length(resamples)){
    cat("Iter ",i,"\n")
    cat("...Prepare train and test sets\n")
    test_id = setdiff(1:nrow(data), resamples[[i]])
    x_train = data[resamples[[i]],]
    y_train = target[resamples[[i]],]
    x_test = data[test_id,] 
    y_test = target[test_id,]
    
    weights = ifelse(y_train == "One", table(y_train)[[2]]/table(y_train)[[1]], 1)
    
    if (any(rownames(x_train) != names(y_train)) | any(rownames(x_test) != names(y_test))){
      message("Samples in train and test sets do not match!\n")
      next
    }
    
    cat("...Train model\n")
    set.seed(seed)
    fit <- caret::train(x = x_train,
                        y = y_train,
                        method="glmnet", 
                        metric=metric,
                        tuneLength = 20,
                        weights = weights,
                        maximize = maximize,
                        trControl=my_control,
                        importance = TRUE)
    
    cat("...Evaluate performance on test set\n")
    pred = predict(fit, x_test, s = "lambda.min", type = "prob")$One
    roc <- roc(response = y_test, predictor = pred, levels = c("Zero","One"))
    auc = auc(roc)
    cat("...Test AUROC is ",auc,"\n")
    perf_list[[i]] = auc
    if (auc > 0.5){
      cat("...Extract important variables\n")
      var_imp = varImp(fit)$importance
      var = var_imp$Overall
      names(var) = rownames(var_imp)
      var = var[order(var, decreasing = T)]
      var = var[var != 0]
      cat("...There are ",length(var)," important variables\n")
      var_list[[i]] = var
    }
    
  }
  
  #Select best features
  var_list <- var_list[lapply(var_list,length)>0]
  if (length(var_list) > 0){
    cat("Aggregate important variable lists\n")
    var_fil = lapply(var_list, function(x) names(x[x>0]))
    set.seed(seed)
    var_rra = aggregateRanks(var_fil)
    var_rra$adjP = var_rra$Score*100
    var_rra$adjP = p.adjust(var_rra$adjP, "fdr")
    var_sel = rownames(var_rra[var_rra$adjP < 0.05,])
    cat("...There are ",length(var_sel)," final variables\n")
  }
  
  return(list(selected_vars = var_sel, performance = perf_list, var_list = var_list))
  
}

genomics_create_samples <- function(
    # creates input files for the genomics feature selection in the working directory
  fam_file = NULL, # path to the .fam file
  fs_samples = NULL, # a vector of samples for Genomics Feature Selection 
  n_resamplings = 100,
  p = 0.8,
  seed = 123
){
  set.seed(seed)
  
  fam = read.csv(fam_file, sep = ' ')
  fam = fam[fam[[2]] %in% fs_samples, ]
  
  
  # making sample files
  samples = createDataPartition(fam[[2]], times=n_resamplings, p=p)
  dir.create('genomics_samples')
  for (i in 1:n_resamplings) {
    fid = fam[samples[[i]], ][[1]]
    iid = fam[samples[[i]], ][[2]]
    write.table(data.frame(V1 = fid, V2 = iid),
                file = file.path(getwd(), paste0("genomics_samples/samples_",i,".txt")),
                col.names = F, row.names = F, quote = F, sep = "\t")
  }
  fid = fam[[1]]
  iid = fam[[2]]
  write.table(data.frame(V1 = fid, V2 = iid),
              file = file.path(getwd(), "genomics_samples/samples_fs.txt"),
              col.names = F, row.names = F, quote = F, sep = "\t")
}

methylomics_create_samples <- function(
    # creates input files for the methylomics feature selection in the working directory
  fs_samples = NULL, # a vector of samples for Methylomics Feature Selection 
  n_resamplings = 100,
  p = 0.8,
  seed = 123
){
  set.seed(seed)
  
  # making sample files
  samples = createDataPartition(fs_samples, times=n_resamplings, p=p)
  dir.create('methylomics_samples')
  saveRDS(fs_samples, file=file.path(getwd(), "methylomics_samples/samples_fs.rds"
  ))
  for (i in 1:n_resamplings) {
    saveRDS(fs_samples[samples[[i]]], file=file.path(getwd(), paste0("methylomics_samples/samples_", i, ".rds")))
  }
}
