#' Feature Selection Methylomics Data
#'
#' @description Feature selection for methylomics data using GSEA.
#'
#' @param trainIDs Set of training IDs from the `data_partitioning()` function
#'   for the feature selection.
#' @param testIDs Set of testing IDs from the `data_partitioning()` function for
#'   the feature selection.
#' @param dataIDs A data.frame with samples as rows and the data modalities as
#'   columns. It holds the data IDs of a sample for each modality. If for a
#'   sample there is no data for a modality, it has to be indicated by NA.
#' @param phenotypeIDs A data.frame with samples as rows and sample IDs as row
#'   names. Columns are phenotypes of interest. If for a sample no information
#'   about a phenotype is available, it has to be indicated by NA.
#' @param methylomicsData A table holding the data for the methylomics.
#' @param lowestLevelPathways A list of pathways used in GSEA.
#' @param seed The seed used for random number generation. Using the same seed
#'   ensures reproducibility.
#'
#' @return Returns methylomicsData subset to the selected features.
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
