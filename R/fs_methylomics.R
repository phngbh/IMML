#' Feature Selection Methylomics Data
#'
#' @description The function for the feature selection for methylomics data.
#'
#' @param train_IDs Set of training IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param test_IDs Set of testing IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param data_IDs A table holding all the clinical IDs with the respective IDs
#' for each modality.
#' @param phenotype_IDs A table holding the clinical IDs with a variable indicating,
#' if the disease occurred or not.
#' @param methylomics_data A table holding the data for the methylomics.
#' @param low_lev_path The lowest level pathways.
#' @param seed The possibility to change the seed for the function.
#'
#' @return Returns the methylomics data table with the selected features.
#'
#' @author Ulrich Asemann

fs_methylomics <- function(train_IDs,
                           test_IDs,
                           data_IDs,
                           phenotype_IDs,
                           methylomics_data,
                           low_lev_path,
                           seed = 123) {
  # Selecting the IDs
  train_methylomics_IDs <-
    train_IDs$`Training Feature Selection IDs`$Methylomics
  test_methylomics_IDs <-
    test_IDs$`Testing Feature Selection IDs`$Methylomics

  # Data frame of phenotypes with all used IDs and their inc3 value
  samples <- unlist(train_methylomics_IDs) %>% unique()
  meth_IDs <-
    data_IDs[match(samples, data_IDs$Clinical), ] %>% dplyr::select("Methylomics") %>% unlist()

  info <-
    phenotype_IDs[as.character(samples), "inc3", drop = FALSE]

  rownames(info) <- meth_IDs

  # Create a model matrix
  info <- info %>% transmute(inc3 = as.character(inc3))

  modmatrix <- model.matrix(~ 0 + ., data = info)

  # Select data, fit, contrast
  data <- methylomics_data[, rownames(modmatrix)]

  fit <- lmFit(data, modmatrix)

  contrast <-
    makeContrasts(inc31 - inc30, levels = colnames(coef(fit)))

  # Temporary data
  tmp <- contrasts.fit(fit, contrast)
  tmp <- eBayes(tmp)

  # topTable function
  topde <-
    topTable(tmp, sort.by = "P", n = Inf)

  # dualmap450kEID in package champ
  # data still missing, check packages from phong for file
  data("dualmap450kEID")

  # geneset_reactome
  geneset_reactome <- reactomePathways(names(mapEIDto450k.lv))
  geneset_reactome <-
    geneset_reactome[intersect(names(geneset_reactome), low_lev_path)]

  # for-loop
  for (p in names(geneset_reactome)) {
    probes <-
      mapEIDto450k.lv[intersect(geneset_reactome[[p]], names(mapEIDto450k.lv))] %>% unlist() %>% unique()
    na_genes <- setdiff(geneset_reactome[[p]], names(mapEIDto450k.lv))
    if (length(probes) > 0) {
      geneset_reactome[[p]] = c(probes, na_genes)
    } else {
      geneset_reactome[[p]] = NULL
    }
  }

  # ranklist
  ranklist <- topde$t
  names(ranklist) <- rownames(topde)
  ranklist <- sort(ranklist)
  set.seed(993)

  # GSEA
  fgseaRes <- fgsea(
    pathways = geneset_reactome,
    stats    = ranklist,
    minSize  = 15,
    maxSize  = 200
  ) %>% arrange(pval) %>% filter(padj < 0.1)

  fgseaRes <-
    fgseaRes %>% mutate(
      leadingEdge = mapIdsList(
        x = org.Hs.eg.db,
        keys = leadingEdge,
        keytype = "ENTREZID",
        column = "SYMBOL"
      )
    )

  # Save data
  saveRDS(fgseaRes, "meth_gsea_final.rds")

  edge_mnsi3 <- fgseaRes$leadingEdge %>% unlist() %>% unique()
  data_selected <- dat[edge_mnsi3,]
  saveRDS(data_selected,"methylomics_selected.rds")

  return("Feature selection Methylomics done!")
}
