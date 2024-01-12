#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(RobustRankAggreg))
suppressMessages(library(fgsea))

args = commandArgs(trailingOnly = TRUE)

setwd(as.character(args[[2]]))

read_file_to_named_list <- function(file_path) {
  # Read the lines of the file
  lines <- readLines(file_path)
  
  # Initialize an empty list to store the data
  named_list <- list()
  
  # Loop through each line to process the data
  for (line in lines) {
    # Remove leading and trailing whitespace
    line <- trimws(line)
    
    # Split the line by space to get the elements
    elements <- strsplit(line, " ")[[1]]
    
    # The first element is the name, and the rest are the IDs
    set_name <- elements[1]
    ids <- elements[-1]
    
    # Add the IDs as character vector to the list with set name as the key
    named_list[[set_name]] <- as.character(ids)
  }
  
  return(named_list)
}

cat("Extract significant pathways\n")
pwlist <- list()
for (i in 1:args[[3]]) {
  res <-
    read.table(
      file.path(
        getwd(),
        "magma_geneset",
        paste0(as.character(args[[1]]), "_", i, ".gsa.out")
      ),
      header = T,
      sep = "",
      skip = 4
    )
  res$adjP <- p.adjust(res$P, method = "fdr")
  res <- arrange(res, adjP) %>% filter(adjP < 0.1)
  pwlist[[i]] <- as.character(res$FULL_NAME)
}

pwlist <- pwlist[lapply(pwlist, length) > 0]
if (length(pwlist) > 3) {
  cat("Aggregate pathways into a single list\n")
  set.seed(993)
  pwlist_agg = aggregateRanks(pwlist)
  pwlist_agg$adjP = pwlist_agg$Score * length(pwlist)
  pwlist_agg$adjP = p.adjust(pwlist_agg$adjP, method = "fdr")
  toppw = filter(pwlist_agg, adjP < 0.05)$Name %>% as.character()
} else {
  cat("There are less than three significant pathways to aggregate\n")
}

if (length(toppw) > 0) {
  write.table(
    data.frame(V1 = toppw),
    file.path(getwd(), paste0(
      as.character(args[[1]]), "top_pathways.txt"
    )),
    col.names = F,
    row.names = F,
    quote = F
  )
  cat("Do final GSEA on most significant pathways\n")
  gene_df <-
    read.table(file.path(getwd(), "magma_gene", paste0(as.character(args[[1]]), ".genes.out")),
               header = T,
               sep = "") %>%
    arrange(ZSTAT)
  ranklist <- gene_df$ZSTAT
  names(ranklist) <- gene_df$GENE
  geneset <- read_file_to_named_list(args[[4]])
  geneset <- geneset[toppw]
  set.seed(993)
  fgseaRes <- fgsea(
    pathways = geneset,
    stats    = ranklist,
    #minSize  = 15,
    maxSize  = 200,
    nperm = 5000
  ) %>% arrange(pval)
  edge_genes <- fgseaRes$leadingEdge %>% unlist() %>% unique()
  
  if (length(edge_genes) > 0) {
    cat("Save leading edge genes and their associated SNPs\n")
    write.table(
      data.frame(V1 = edge_genes),
      file.path(getwd(), paste0(
        as.character(args[[1]]), "_leadingEdge_genes.txt"
      )),
      col.names = F,
      row.names = F,
      quote = F
    )
    gene_snps <-
      read.table(file.path(getwd(), 'annotates', "gene_annot.txt")) %>% filter(V1 %in% edge_genes)
    
    write.table(
      data.frame(V1 = unique(gene_snps$V2)),
      file.path(getwd(), paste0(
        as.character(args[[1]]), "_leadingEdge_SNPs.txt"
      )),
      col.names = F,
      row.names = F,
      quote = F
    )
  } else {
    cat("There are no leading edge genes\n")
  }
} else {
  cat("There are no significant pathways after aggregation\n")
}
