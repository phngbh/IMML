#' Title
#'
#' @param lowestLevelPathways
#' @param genomicsAnnotation
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
FsGenomics <- function(lowestLevelPathways,
                       genomicsAnnotation,
                       seed = 123) {
  # # Geneset reactome
  # genesetReactome <-
  #   reactomePathways(levels(as.factor(as.character(
  #     genomicsAnnotation$V1
  #   ))))
  # genesetReactome <-
  #   genesetReactome[intersect(names(genesetReactome), lowestLevelPathways)]
  #
  # Loop through all possible pathways
  # length(names(genesetReactome)) = 1990
  #
  # for (p in names(genesetReactome)) {
  #   iter <- match(p, names(genesetReactome))
  #   cat("Pathway:", iter, "\n")
  #
  #   gen <-
  #     intersect(levels(as.factor(as.character(
  #       genomicsAnnotation$V1
  #     ))), genesetReactome[[p]])
  #   df <- filter(genomicsAnnotation, V1 %in% as.integer(gen))
  #   na_genes <- setdiff(genesetReactome[[p]], levels(as.factor(as.character(geno_annot$V1))))
  #   if (length(gen) > 0) {
  #     genesetReactome[[p]] <- c(df$V2 %>% unique(), na_genes)
  #   } else {
  #     genesetReactome[[p]] <- NULL
  #     cat("no pathway")
  #   }
  # }
  #
  # return(genesetReactome)
  # # Data was saved in folder "data"

  # Load data
  data("genomics_geneset_reactome")
  genesetReactome <- genomics_geneset_reactome

  # Define variables
  gsea = list()

  for (i in 1:100) {
    cat("Iter", i, "\n")

    # Loading the data for the iteration
    gwasTmp <-
      read.csv(
        paste0(
          "C:/Users/ulric/Desktop/Arbeit/data/IMLdata/download/fs_genomics/samples_",
          i,
          ".inc3.assoc.logistic"
        ),
        sep = ""
      ) %>%
      filter(SNP %in% genomicsAnnotation$V2)

    # Select the needed data
    gwasTmp$GENE <-
      genomicsAnnotation$V1[match(gwasTmp$SNP, genomicsAnnotation$V2)]
    gwasTmpUnq <- gwasTmp %>% distinct(STAT, GENE, .keep_all = T)
    ranklist <- gwasTmpUnq$STAT
    names(ranklist) <- gwasTmpUnq$SNP
    ranklist <- sort(ranklist)

    # Set seed
    set.seed(seed)

    # Do fgsea
    fgseaResTmp <- fgsea(
      pathways = genesetReactome,
      stats    = ranklist,
      eps = 0,
      minSize = 15,
      maxSize  = 500
    ) %>% arrange(pval) #%>% filter(padj < 5e-8)

    if (nrow(fgseaResTmp) == 0) {
      message("No significant pathway found\n")
      next
    }

    fgseaResTmp$leadingEdge_gen = lapply(fgseaResTmp$leadingEdge,
                                          function(x)
                                            return(gwasTmpUnq$GENE[match(x, gwasTmpUnq$SNP)] %>%
                                                     unique()))
    gsea[[i]] = fgseaResTmp
  }

  # Pathways
  pathways = list()
  for (i in 1:length(gsea)) {
    if (is.null(gsea[[i]])) {
      next
    }
    pathways[[i]] = gsea[[i]]$pathway
  }

  pwlistAgg <- aggregateRanks(pathways)
  pwlistAgg$adjP <- pwlistAgg$Score * length(pathways)
  pwlistAgg$adjP <- p.adjust(pwlistAgg$adjP, method = "fdr")
  toppw <- rownames(pwlistAgg[pwlistAgg$adjP < 0.05,])
  # if adjP doesn't work, use code below
  # toppw <- rownames(pwlistAgg[pwlistAgg$adjP,])

  # Read data
  gwas = read.csv(
    paste0(
      "C:/Users/ulric/Desktop/Arbeit/data/IMLdata/download/fs_genomics/all.inc3.assoc.logistic"
    ),
    sep = ""
  ) %>%
    filter(SNP %in% genomicsAnnotation$V2)

  # Do gwas
  gwas$GENE <-
    genomicsAnnotation$V1[match(gwas$SNP, genomicsAnnotation$V2)]
  gwasMnsi3Unq <- gwas %>% distinct(STAT, GENE, .keep_all = T)
  ranklist <- gwasMnsi3Unq$STAT
  names(ranklist) <- gwasMnsi3Unq$SNP
  ranklist <- sort(ranklist)

  # Final fgsea
  fgseaRes <- fgsea(
    pathways = genesetReactome[toppw],
    stats    = ranklist,
    eps = 0,
    minSize = 15,
    maxSize  = 500
  ) %>% arrange(pval)

  fgseaRes$leadingEdge_gen <- lapply(fgseaRes$leadingEdge,
                                     function(x)
                                       return(gwasMnsi3Unq$GENE[match(x, gwasMnsi3Unq$SNP)] %>%
                                                unique()))

  fgseaRes$leadingEdge_s <- lapply(fgseaRes$leadingEdge_gen,
                                   function(x)
                                     sort(x) %>% paste0(collapse = " ")) %>% unlist()
  fgseaRes$similarity <- sapply(fgseaRes$leadingEdge_s,
                                function(x)
                                  agrep(x, fgseaRes$leadingEdge_s, max = 0.2))
  fgseaResFil <- distinct(fgseaRes, similarity, .keep_all = T)

  edge <- fgseaRes$leadingEdge %>% unlist() %>% unique()

  # Write data into file
  write.table(
    data.frame(V1 = edge),
    file = "edge_snp.txt",
    col.names = F,
    row.names = F,
    quote = F
  )


  return("Feature selection genomics done!")

}
