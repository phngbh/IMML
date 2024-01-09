#!/usr/bin/env Rscript

library(dplyr)
library(RobustRankAggreg)
library(fgsea)
library(reactome.db)

args = commandArgs(trailingOnly = TRUE)

setwd(as.character(args[[1]]))

cat("Extract significant genes\n")
sig_gene_list <- list()
for (i in 1:args[[3]]) {
  file <-
    file.path(getwd(), 'magma_gene', paste0(args[[2]], "_", i, ".genes.out"))
  if (file.exists(file)) {
    gene_stat <- read.table(file, header = T)
    gene_stat$ADJP <- p.adjust(gene_stat$P, method = "fdr")
    gene_stat <- arrange(gene_stat, ADJP)
    sig_genes <- gene_stat$GENE[gene_stat$ADJP < 0.05]
    sig_gene_list[[i]] <- sig_genes
  } else {
    message("...File not found for iteration ", i, "\n")
  }
}

sig_gene_list <- sig_gene_list[lapply(sig_gene_list, length) > 0]

cat("Aggregate lists of significant genes into a single list\n")
if (length(sig_gene_list) <= 3) {
  cat("...There are less than 3 significant gene lists => no aggregate\n")
  
} else {
  set.seed(993)
  print(sig_gene_list)
  gene_agg = aggregateRanks(glist = sig_gene_list)
  gene_agg$adjP = gene_agg$Score * length(sig_gene_list)
  gene_agg$adjP = p.adjust(gene_agg$adjP, method = "fdr")
  aggregated_genes = rownames(filter(gene_agg, adjP < 0.05))
  cat(paste0(
    "...There are ",
    length(aggregated_genes),
    " aggregated genes\n"
  ))
  if (length(aggregated_genes) > 0) {
    saveRDS(aggregated_genes,
            paste0(args[[2]], "_significant_genes.rds"))
    cat("...Extract associated SNPs\n")
    gene_annot <-
      read.table("gene_annot.txt") %>% dplyr::filter(V1 %in% aggregated_genes)
    write.table(
      data.frame(V1 = unique(gene_annot$V2)),
      file.path(getwd(), paste0(
        as.character(args[[2]]), "_associated_SNPs.txt"
      )),
      col.names = F,
      row.names = F,
      quote = F
    )
  }
}
