fs_genomics <- function(low_lev_path,
                        genomics_annot,
                        seed = 123) {
  # Geneset reactome
  geneset_reactome <-
    reactomePathways(levels(as.factor(as.character(
      genomics_annot$V1
    ))))
  geneset_reactome <-
    geneset_reactome[intersect(names(geneset_reactome), low_lev_path)]

  # return(length(names(geneset_reactome)))


  # names(geneset_reactome)
  for (p in 1:5) {
    gen <-
      intersect(levels(as.factor(as.character(
        genomics_annot$V1
      ))), geneset_reactome[[p]])
    df = filter(genomics_annot, V1 %in% as.integer(gen))
    na_genes = setdiff(geneset_reactome[[p]], levels(as.factor(as.character(geno_annot$V1))))
    if (length(gen) > 0) {
      geneset_reactome[[p]] = c(df$V2 %>% unique(), na_genes)
    } else {
      geneset_reactome[[p]] = NULL
    }
  }

  # return(geneset_reactome)


  # Building gsea list, code line 18

  gsea = list()
  # i in 1:100
  for (i in 1:1) {
    cat("Iter", i, "\n")
    gwas_tmp <-
      read.csv(
        paste0(
          "C:/Users/ulric/Desktop/Arbeit/data/IMLdata/download/fs_genomics/samples_",
          i,
          ".inc3.assoc.logistic"
        ),
        sep = ""
      ) %>%
      filter(SNP %in% genomics_annot$V2)

    gwas_tmp$GENE <-
      genomics_annot$V1[match(gwas_tmp$SNP, genomics_annot$V2)]
    gwas_tmp_unq <- gwas_tmp %>% distinct(STAT, GENE, .keep_all = T)
    ranklist <- gwas_tmp_unq$STAT
    names(ranklist) <- gwas_tmp_unq$SNP
    ranklist <- sort(ranklist)

    # if (i == 2){
    #   return(gwas_tmp)
    # }

    set.seed(seed)

    fgseaRes_tmp <- fgsea(
      pathways = geneset_reactome,
      stats    = ranklist,
      eps = 0,
      # minSize = 15
      minSize  = 1,
      maxSize  = 500
    ) %>% arrange(pval) #%>% filter(padj < 5e-8)

    if (nrow(fgseaRes_tmp) == 0) {
      message("No significant pathway found\n")
      next
    }

    fgseaRes_tmp$leadingEdge_gen = lapply(fgseaRes_tmp$leadingEdge,
                                          function(x)
                                            return(gwas_tmp_unq$GENE[match(x, gwas_tmp_unq$SNP)] %>%
                                                     unique()))
    gsea[[i]] = fgseaRes_tmp

    # return(fgseaRes_tmp)

  }

  # Pathway list
  pathways = list()
  for (i in 1:length(gsea)) {
    if (is.null(gsea[[i]])) {
      next
    }
    pathways[[i]] = gsea[[i]]$pathway
  }

  # return(pathways)

  pwlist_agg <- aggregateRanks(pathways)
  pwlist_agg$adjP <- pwlist_agg$Score * length(pathways)
  pwlist_agg$adjP <- p.adjust(pwlist_agg$adjP, method = "fdr")
  # toppw <- rownames(pwlist_agg[pwlist_agg$adjP < 0.05,])
  toppw <- rownames(pwlist_agg[pwlist_agg$adjP,])

  gwas = read.csv(
    paste0(
      "C:/Users/ulric/Desktop/Arbeit/data/IMLdata/download/fs_genomics/all.inc3.assoc.logistic"
    ),
    sep = ""
  ) %>%
    filter(SNP %in% genomics_annot$V2)

  # return(gwas)

  gwas$GENE <- genomics_annot$V1[match(gwas$SNP, genomics_annot$V2)]
  gwas_unq <- gwas %>% distinct(STAT, GENE, .keep_all = T)
  ranklist <- gwas_unq$STAT
  names(ranklist) <- gwas_unq$SNP
  ranklist <- sort(ranklist)

  # return(toppw)
  # return(geneset_reactome)

  # Final gsea
  fgseaRes <- fgsea(
    pathways = geneset_reactome[toppw[1]],
    stats    = ranklist,
    eps = 0,
    # minSize = 15
    minSize  = 1,
    maxSize  = 500
  ) %>% arrange(pval)

  return(fgseaRes)
  fgseaRes$leadingEdge_gen <- lapply(fgseaRes$leadingEdge,
                                    function(x)
                                      return(gwas_mnsi3_unq$GENE[match(x, gwas_mnsi3_unq$SNP)] %>%
                                               unique()))

  fgseaRes$leadingEdge_s <- lapply(fgseaRes$leadingEdge_gen,
                                  function(x)
                                    sort(x) %>% paste0(collapse = " ")) %>% unlist()
  fgseaRes$similarity <- sapply(fgseaRes$leadingEdge_s,
                               function(x)
                                 agrep(x, fgseaRes$leadingEdge_s, max = 0.2))
  fgseaRes_fil <- distinct(fgseaRes, similarity, .keep_all = T)

  edge <- fgseaRes$leadingEdge %>% unlist() %>% unique()
  write.table(
    data.frame(V1 = edge),
    file = "edge_snp.txt",
    col.names = F,
    row.names = F,
    quote = F
  )


  return("Feature selection genomics done!")

}
