#' Feature Selection Proteomics Data
#'
#' @description The function for the feature selection for proteomics data.
#'
#' @param trainIDs Set of training IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param testIDs Set of testing IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param dataIDs A table holding all the clinical IDs with the respective IDs
#' for each modality.
#' @param phenotypeIDs A table holding the clinical IDs with a variable indicating,
#' if the disease occurred or not.
#' @param lowestLevelPathways The lowest level pathways.
#' @param somamaerInfoEdited An annotation file with additional information
#' for the proteomics.
#' @param seed The possibility to change the seed for the function.
#' @param proteomicsData A table holding the data for the proteomics.
#'
#' @return Returns the proteomics data table with the selected features.
#'
#' @author Ulrich Asemann
#'
fs_proteomics <- function(trainIDs,
                          testIDs,
                          dataIDs,
                          phenotypeIDs,
                          proteomicsData,
                          lowestLevelPathways,
                          somamaerInfoEdited,
                          seed = 123) {
  # Selecting the IDs
  trainProteomicsIDs <-
    trainIDs$`Training Feature Selection IDs`$Proteomics
  testProteomicsIDs <-
    testIDs$`Testing Feature Selection IDs`$Proteomics

  # Data frame of phenotypes with all used IDs and their inc3 value
  samples <- unlist(trainProteomicsIDs) %>% unique()
  protIDs <-
    dataIDs[match(samples, dataIDs$Clinical),] %>% dplyr::select("Proteomics") %>% unlist()

  info <-
    phenotypeIDs[as.character(samples), "inc3", drop = FALSE]
  rownames(info) <- protIDs

  # Creating a model matrix
  info <- info %>% transmute(inc3 = as.character(inc3))

  modmatrix <- model.matrix( ~ 0 + ., data = info)

  # Select data, fit, contrast
  data <- proteomicsData[, rownames(modmatrix)]

  fit <- lmFit(data, modmatrix)

  contrast <-
    makeContrasts(inc31 - inc30, levels = colnames(coef(fit)))

  # Temporary data
  tmp <- contrasts.fit(fit, contrast)
  tmp <- eBayes(tmp)

  # topTable function
  topde <-
    topTable(tmp, sort.by = "P", n = Inf)

  topde <-
    topde %>%
    dplyr::mutate(Gene = rownames(topde)) %>%
    dplyr::mutate(Name = somamaerInfoEdited$EntrezGeneSymbol
                  [match(.$Gene, somamaerInfoEdited$SeqId)],
                  EntrezID = somamaerInfoEdited$EntrezGeneID
                  [match(.$Gene, somamaerInfoEdited$SeqId)])

  # Select from topde
  topdeTmp <- dplyr::select(topde, EntrezID, P.Value) %>% na.omit()
  tmp <- tapply(topdeTmp$P.Value, topdeTmp$EntrezID, min)

  # Create lists
  ranklist <- vector(mode = "numeric", length = length(tmp))
  names(ranklist) <- names(tmp)

  ilmnlist <- vector(mode = "character", length = length(tmp))
  names(ilmnlist) <- names(tmp)

  # Add values to the lists
  tmp2 <-
    dplyr::filter(topde, EntrezID %in% names(tmp) &
                    P.Value %in% tmp)

  tmp2 <- tmp2 %>% dplyr::select(EntrezID, t, Gene)
  ranklist[as.character(tmp2$EntrezID)] <- tmp2$t
  ilmnlist[as.character(tmp2$EntrezID)] <- tmp2$Gene


  # Build dataframe with the lists
  genelist = data.frame(Entrez = names(ranklist),
                        Probe = ilmnlist,
                        t = ranklist)

  # Get the pathways
  ranklist = sort(ranklist)
  genesetReactome = reactomePathways(names(ranklist))

  # Add this part for the pathways,
  # when 100 iterations are implemented in the code
  # genesetReactome = genesetReactome[intersect(names(genesetReactome), toppw)]

  # Set seed
  set.seed(993)

  # Fast GSEA
  fgseaRes <- fgsea(
    pathways = genesetReactome,
    stats    = ranklist,
    minSize  = 5,
    maxSize  = 200
  ) %>% arrange(pval) %>% filter(padj < 0.1)

  fgseaRes <-
    fgseaRes %>% mutate(
      leadingEdge = mapIdsList(
        x = org.Hs.eg.db,
        keys = fgseaRes$leadingEdge,
        keytype = "ENTREZID",
        column = "SYMBOL"
      )
    )

  # Save results from fgseaRes
  saveRDS(fgseaRes, "prot_gsea_final.rds")

  # Extract selected features for training
  edge <- fgseaRes$leadingEdge %>% unlist() %>% unique()
  edgeEntrez <- AnnotationDbi::select(org.Hs.eg.db,
                                       keys = edge,
                                       columns = "ENTREZID",
                                       keytype = "SYMBOL")
  probe <- genelist$Probe[genelist$Entrez %in% edgeEntrez$ENTREZID]
  dataSelected <- proteomicsData[probe, , drop = F]

  # Save the selected data
  saveRDS(dataSelected, "proteomics_selected.rds")

  return("Feature selection proteomics done!")
}
