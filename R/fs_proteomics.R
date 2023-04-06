#' Title
#'
#' @param train_IDs
#' @param test_IDs
#' @param data_IDs
#' @param phenotype_IDs
#' @param low_lev_path
#' @param somamaer_info_edited
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
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
    train_IDs$`Feature Selection IDs`$Proteomics
  test_proteomics_IDs <-
    test_IDs$`Feature Selection IDs`$Proteomics

  # Data frame of phenotypes with all used IDs and their inc3 value
  samples <- unlist(train_proteomics_IDs) %>% unique()
  prot_IDs <-
    data_IDs[match(samples, data_IDs$Clinical), ] %>% dplyr::select("Proteomics") %>% unlist()

  info <-
    phenotype_IDs[as.character(samples), "inc3", drop = FALSE]
  rownames(info) <- prot_IDs


  # Creating a model matrix
  info <- info %>% transmute(inc3 = as.character(inc3))

  modmatrix <- model.matrix(~ 0 + ., data = info)

  # Select data, fit, contrast
  data <- proteomics_data[, rownames(modmatrix)]

  fit <- lmFit(dat, modmatrix)

  contrast <-
    makeContrasts(inc31 - inc30, levels = colnames(coef(fit)))

  # Temporary data
  tmp <- contrasts.fit(fit, contrast)
  tmp <- eBayes(tmp)

  # topTable function
  topde <-
    topTable(tmp, sort.by = "P", n = Inf) %>%
    mutate(Name = somamaer_info_edited$EntrezGeneSymbol
           [match(.$Gene, rownames(somamaer_info_edited))],
           EntrezID = somamaer_info_edited$EntrezGeneID
           [match(.$Gene, rownames(somamaer_info_edited))])

  topde.tmp <- dplyr::select(topde, EntrezID, P.Value) %>% na.omit()
  tmp <- tapply(topde.tmp$P.Value, topde.tmp$EntrezID, min)

  # Build lists
  ranklist <- vector(mode = "numeric", length = length(tmp))
  names(ranklist) <- names(tmp)

  ilmnlist <- vector(mode = "character", length = length(tmp))

  for (i in length(tmp)) {
    # Get values
    t <-
      dplyr::filter(topde, EntrezID == names(tmp)[i] &
                      P.Value == tmp[i]$t)
    il <-
      dplyr::filter(topde, EntrezID == names(tmp)[i] &
                      P.Value == tmp[i]$Gene)

    # Add to list
    ranklist[i] <- t
    ilmnlist[i] <- il
  }

  # Build dataframe
  genelist = data.frame(Entrez = names(ranklist),
                        Probe = ilmnlist,
                        t = ranklist)
  ranklist = sort(ranklist)
  geneset_reactome = reactomePathways(names(ranklist))
  geneset_reactome = geneset_reactome[intersect(names(geneset_reactome), toppw)]
  set.seed(993)

  # GSEA
  fgseaRes <- fgsea(
    pathways = geneset_reactome,
    stats    = ranklist,
    minSize  = 5,
    maxSize  = 200
  ) %>% arrange(pval) %>% filter(padj < 0.1)
  fgseaRes[, leadingEdge := mapIdsList(
    x = org.Hs.eg.db,
    keys = leadingEdge,
    keytype = "ENTREZID",
    column = "SYMBOL"
  )]

  # Save results
  saveRDS(fgseaRes, "prot_gsea_final.rds")

  # Extract selected features for training
  edge <- fgseaRes$leadingEdge %>% unlist() %>% unique()
  edge_entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                       keys = edge,
                                       columns = "ENTREZID",
                                       keytype = "Symbol")
  probe <- genelist$Probe[genelist$Entrez %in% edge_entrez$ENTREZID]
  dat_selected <- proteomics_data[probe, , drop = F]

  # Save selected data
  saveRDS(dat_selected, "proteomics_selected.rds")



  return("Feature selection proteomics done!")

}
