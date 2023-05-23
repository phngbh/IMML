#' Feature Selection Proteomics Data
#'
#' @description The function for the feature selection for proteomics data.
#'
#' @param train_IDs Set of training IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param test_IDs Set of testing IDs from the `data_partitioning()` function
#' for the feature selection.
#' @param data_IDs A table holding all the clinical IDs with the respective IDs
#' for each modality.
#' @param phenotype_IDs A table holding the clinical IDs with a variable indicating,
#' if the disease occurred or not.
#' @param low_lev_path The lowest level pathways.
#' @param somamaer_info_edited An annotation file with additional information
#' for the proteomics.
#' @param seed The possibility to change the seed for the function.
#' @param proteomics_data A table holding the data for the proteomics.
#'
#' @return Returns the proteomics data table with the selected features.
#'
#' @author Ulrich Asemann
#'
fs_proteomics <- function(train_IDs,
                          test_IDs,
                          data_IDs,
                          phenotype_IDs,
                          proteomics_data,
                          low_lev_path,
                          somamaer_info_edited,
                          seed = 123) {
  # Selecting the IDs
  train_proteomics_IDs <-
    train_IDs$`Training Feature Selection IDs`$Proteomics
  test_proteomics_IDs <-
    test_IDs$`Testing Feature Selection IDs`$Proteomics

  # Data frame of phenotypes with all used IDs and their inc3 value
  samples <- unlist(train_proteomics_IDs) %>% unique()
  prot_IDs <-
    data_IDs[match(samples, data_IDs$Clinical),] %>% dplyr::select("Proteomics") %>% unlist()

  info <-
    phenotype_IDs[as.character(samples), "inc3", drop = FALSE]
  rownames(info) <- prot_IDs

  # Creating a model matrix
  info <- info %>% transmute(inc3 = as.character(inc3))

  modmatrix <- model.matrix( ~ 0 + ., data = info)

  # Select data, fit, contrast
  data <- proteomics_data[, rownames(modmatrix)]

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
    dplyr::mutate(Name = somamaer_info_edited$EntrezGeneSymbol
                  [match(.$Gene, somamaer_info_edited$SeqId)],
                  EntrezID = somamaer_info_edited$EntrezGeneID
                  [match(.$Gene, somamaer_info_edited$SeqId)])

  # Select from topde
  topde_tmp <- dplyr::select(topde, EntrezID, P.Value) %>% na.omit()
  tmp <- tapply(topde_tmp$P.Value, topde_tmp$EntrezID, min)

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
  geneset_reactome = reactomePathways(names(ranklist))

  # Add this part for the pathways,
  # when 100 iterations are implemented in the code
  # geneset_reactome = geneset_reactome[intersect(names(geneset_reactome), toppw)]

  # Set seed
  set.seed(993)

  # Fast GSEA
  fgseaRes <- fgsea(
    pathways = geneset_reactome,
    stats    = ranklist,
    minSize  = 5,
    maxSize  = 200
  ) %>% arrange(pval) #%>% filter(padj < 0.1)

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
  edge_entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                       keys = edge,
                                       columns = "ENTREZID",
                                       keytype = "SYMBOL")
  probe <- genelist$Probe[genelist$Entrez %in% edge_entrez$ENTREZID]
  dat_selected <- proteomics_data[probe, , drop = F]

  # Save the selected data
  saveRDS(dat_selected, "proteomics_selected.rds")

  return("Feature selection proteomics done!")
}
