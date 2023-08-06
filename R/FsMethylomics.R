#' Feature Selection Methylomics Data
#'
#' @description The function for the feature selection for methylomics data.
#'
#' @param trainIDs Set of training IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param testIDs Set of testing IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param dataIDs A table holding all the clinical IDs with the respective IDs
#' for each modality.
#' @param phenotypeIDs A table holding the clinical IDs with a variable indicating,
#' if the disease occurred or not.
#' @param methylomicsData A table holding the data for the methylomics.
#' @param lowestLevelPathways The lowest level pathways.
#' @param seed The possibility to change the seed for the function.
#'
#' @return Returns the methylomics data table with the selected features.
#'
#' @author Ulrich Asemann

FsMethylomics <- function(trainIDs,
                          testIDs,
                          dataIDs,
                          phenotypeIDs,
                          methylomicsData,
                          lowestLevelPathways,
                          seed = 123) {
  # Selecting the IDs
  trainMethylomicsIDs <-
    trainIDs$`Training Feature Selection IDs`$Methylomics
  testMethylomicsIDs <-
    testIDs$`Testing Feature Selection IDs`$Methylomics

  # Data frame of phenotypes with all used IDs and their inc3 value
  samples <- unlist(trainMethylomicsIDs) %>% unique()
  methIDs <-
    dataIDs[match(samples, dataIDs$Clinical), ] %>% dplyr::select("Methylomics") %>% unlist()

  info <-
    phenotypeIDs[as.character(samples), "inc3", drop = FALSE]
  rownames(info) <- methIDs

  # Create a model matrix
  info <- info %>% transmute(inc3 = as.character(inc3))
  modmatrix <- model.matrix(~ 0 + ., data = info)

  # Select data, fit, contrast
  data <- methylomicsData[, rownames(modmatrix)]
  fit <- lmFit(data, modmatrix)
  contrast <-
    makeContrasts(inc31 - inc30, levels = colnames(coef(fit)))

  # Temporary data
  tmp <- contrasts.fit(fit, contrast)
  tmp <- eBayes(tmp)

  # topTable function
  topde <-
    topTable(tmp, sort.by = "P", n = Inf)

  # genesetReactome
  data("dualmap450kEID")
  genesetReactome <- reactomePathways(names(mapEIDto450k.lv))
  genesetReactome <-
    genesetReactome[intersect(names(genesetReactome), lowestLevelPathways)]

  # for-loop
  for (p in names(genesetReactome)) {
    probes <-
      mapEIDto450k.lv[intersect(genesetReactome[[p]], names(mapEIDto450k.lv))] %>% unlist() %>% unique()
    naGenes <-
      setdiff(genesetReactome[[p]], names(mapEIDto450k.lv))
    if (length(probes) > 0) {
      genesetReactome[[p]] = c(probes, naGenes)
    } else {
      genesetReactome[[p]] = NULL
    }
  }

  # ranklist
  ranklist <- topde$t
  names(ranklist) <- rownames(topde)
  ranklist <- ranklist[order(names(ranklist))]

  # seed
  set.seed(993)

  # GSEA
  fgseaRes <- fgsea(
    pathways = genesetReactome,
    stats    = ranklist,
    minSize  = 15,
    maxSize  = 200
  ) %>% arrange(pval) %>% filter(padj < 0.05)

  edgeMnsi3 <- fgseaRes$leadingEdge %>% unlist() %>% unique()
  dataSelected <- methylomicsData[edgeMnsi3,]

  # Save data
  saveRDS(dataSelected, "methylomics_selected.rds")
  saveRDS(fgseaRes, "meth_gsea_final.rds")

  return("Feature selection Methylomics done!")
}
