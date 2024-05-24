#!/usr/bin/env Rscript

## USAGE: Rscript processing_methylation.R <working dir> <phenotype file> <sample folder> <data file> <iterations> <pathway list> <dualmap file> <DM threshold> <target_name>

args = commandArgs(trailingOnly=TRUE)
setwd(as.character(args[[1]]))

suppressMessages(library(dplyr))
suppressMessages(library(haven))
suppressMessages(library(matrixStats))
suppressMessages(library(MatrixGenerics))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(bumphunter))
suppressMessages(library(minfi))
suppressMessages(library(DMRcate))
#suppressMessages(library(ChAMPdata))
suppressMessages(library(ChAMP))
suppressMessages(library(limma))
suppressMessages(library(fgsea))
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
#suppressMessages(library(tidyverse))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(caret))
suppressMessages(library(glmnet))
suppressMessages(library(pROC))
#suppressMessages(library(reactome.db))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(methylGSA))
suppressMessages(library(RobustRankAggreg))

set.seed(993)

cat("Load necessary files and prepare for statistical analysis\n")
cat("...Load phenotype file\n")
# Load phenotype file
phenotype <- readRDS(as.character(args[[2]]))
target <- colnames(phenotype)[1]
cat("......Phenotype file loaded with dimensions: ",dim(phenotype),"\n")

cat("...Load pathways\n")
# Load gene set data
#lowest_level_pathways = readRDS("lowest_level_pathways.rds")
load(as.character(args[[7]]))
#geneset_reactome = reactomePathways(names(mapEIDto450k.lv))
#geneset_reactome <- geneset_reactome[intersect(names(geneset_reactome), lowest_level_pathways)]
pathway_list <- readRDS(as.character(args[[6]]))

cat("...Load data file\n")
# Load processed (filtered and batch corrected) data
beta_processed = readRDS(as.character(args[[4]]))
cat("......Data table loaded with dimensions: ", dim(beta_processed),"\n")

cat("...Load samples\n")
# Load sample file
fs_samples <- readRDS(paste0(args[[3]],"/samples_fs.rds"))
cat("......There are: ", length(fs_samples), " samples loaded\n")

cat("...Extract data and make model matrix for feature selection\n")
# Extract data and model matrix for DMA
int_fs <- intersect(fs_samples, intersect(colnames(beta_processed),rownames(phenotype)))
cat("......There are: ", ncol(beta_processed[,int_fs]), " samples in the feature selection set\n")
phenotype <- phenotype[int_fs,,drop=FALSE]
cat("......Dimensions of filtered phenotype table: ",dim(phenotype),"\n")
phenotype[,target] <- as.factor(phenotype[,target])
f <- reformulate(termlabels = target, intercept = F)
modmatrix <- model.matrix(f, data = phenotype)
# modmatrix <- model.matrix(~ 0 + ., data = phenotype)
cat("......Dimensions of filtered model matrix: ",dim(modmatrix),"\n")

# Initiate result list
res = list(edge = list(), ilmn = list(), roc = list(), auc = list(), pathway = list())

# Get subset number
#i <- as.numeric(args[[5]])

iters = as.numeric(args[[5]])

for (i in 1:iters){
  
  cat("Start feature selection for resample ",i,"\n")
  
  cat("...Extract relevant samples\n")
  subset <- readRDS(paste0(args[[3]],"/samples_",i,".rds"))
  int_sub <- intersect(subset, int_fs)
  modmatrix_sub <- modmatrix[colnames(beta_processed[,int_sub]),,drop = FALSE]
  cat("......Subset has ", length(int_sub), " samples\n")
  
  # Differential methylation analysis
  cat("...DM analysis\n")
  if( any(colnames(beta_processed[,int_sub]) != rownames(modmatrix_sub))){
    stop("......columns of data do not match rows of design matrix => Stop!")
  }
  if( any(table(modmatrix_sub[,1]) < 2) | any(table(modmatrix_sub[,2]) < 2)){
    stop("......One of the factor levels has less than 2 observations => Stop!\n")
  }
  fit = lmFit(log(beta_processed[,int_sub]), modmatrix_sub)
  contrast = makeContrasts(paste0(target,"1-",target,"0"), levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contrast)
  tmp <- eBayes(tmp)
  topdm <- topTable(tmp, sort.by = "P", n = Inf)
  sigdm <- rownames(topdm[topdm$P.Value < as.numeric(args[[8]]),])
  cat("......There are ", length(sigdm), " significant differentially methylated probes\n")
  saveRDS(sigdm, file.path(getwd(), paste0("sigDM_",target,"_",i,".rds")))
  ranklist = topdm$P.Value
  names(ranklist) = rownames(topdm)
  ranklist = sort(ranklist)
  
  # Gene set enrichment analysis
  cat("...GSEA\n")
  gseaRes <- methylRRA(cpg.pval = ranklist, method = "GSEA", GS.list = pathway_list, GS.idtype = "ENTREZID", minsize = 3, maxsize = 200) %>%
    dplyr::filter(padj < 0.05)
  if (nrow(gseaRes) == 0){
    cat("......There is no significantly enriched pathway => Stop!\n")
    res[["edge"]] = NULL
    res[["ilmn"]] = NULL
    res[["pathway"]] = NULL
  } else {
    edge_tmp <- gseaRes$core_enrichment %>% paste(collapse = "/") %>% strsplit(split = "/") %>% unlist() %>% unique()
    eid <- mapIdsList(x=org.Hs.eg.db, keys=edge_tmp,keytype="SYMBOL", column="ENTREZID") %>% unlist() %>% unname()
    ilmn <- mapEIDto450k.lv[eid] %>% unlist() %>% unique()
    ilmn <- intersect(ilmn, rownames(beta_processed))
    res[["edge"]] = edge_tmp
    res[["ilmn"]] = ilmn
    res[["pathway"]] = gseaRes$ID
    
    cat("...Lasso\n")
    probelist_tmp = res$ilmn
    oob = setdiff(colnames(beta_processed[,int_fs]), colnames(beta_processed[,int_sub]))
    x_train = beta_processed[,int_sub][probelist_tmp,] %>% t()
    y_train = phenotype[colnames(beta_processed[,int_sub]),target] %>% droplevels()
    up = upSample(x = x_train, y = y_train)
    x_train_up = up[,-ncol(up)]
    y_train_up = up$Class
    x_test = beta_processed[,int_fs][probelist_tmp,oob] %>% t()
    y_test = phenotype[oob,target] %>% droplevels()
    if (nlevels(y_train_up) < 0 | nlevels(y_test) < 2) {
      gsea_inc1$roc = NA
      gsea_inc1$auc = NA
      cat("......The labels in train or test set has only one level => Stop!\n")
    } else {
      fit_tmp = cv.glmnet(x = as.matrix(x_train_up), y = y_train_up, alpha = 1, family ="binomial",
                          nfolds = 3, type.measure = "auc")
      pred_prob = predict(fit_tmp, x_test, s = "lambda.min", type = "response")
      roc <- roc(response = y_test, predictor = pred_prob[,1], levels = c("0","1"))
      auc = auc(roc)
      res[["roc"]] = roc
      res[["auc"]] = auc
    }
  }
  
  saveRDS(res, paste0("fsRes_",target,"_",i,".rds"))
}

tmp_gsea_genes <- list()
tmp_dm_probes <- list()
for (i in 1:iters){
  .gseaRes <- readRDS(paste0("fsRes_",target,"_",i,".rds"))
  if (length(.gseaRes$auc) != 0 && .gseaRes$auc > 0.5){
    tmp_gsea_genes[[i]] <- .gseaRes$edge
  }
  tmp_dm_probes[[i]] <- readRDS(paste0("sigDM_",target,"_",i,".rds"))
}

tmp_gsea_genes <- tmp_gsea_genes[lapply(tmp_gsea_genes, length) > 0]
tmp_genes_agg = aggregateRanks(tmp_gsea_genes)
tmp_genes_agg$adjP = tmp_genes_agg$Score*length(tmp_gsea_genes)
tmp_genes_agg$adjP = p.adjust(tmp_genes_agg$adjP, method = "fdr")
aggregated_genes = rownames(filter(tmp_genes_agg, adjP < 0.05))

if (length(aggregated_genes) == 0){
  tmp_associated_probes = character(0)
} else {
  tmp <- mapIds(org.Hs.eg.db, aggregated_genes, 'ENTREZID', 'SYMBOL')
  tmp_associated_probes <- mapEIDto450k.lv[as.character(unique(tmp))] %>% unlist() %>% unique()
}

tmp_probes_agg = aggregateRanks(tmp_dm_probes)
tmp_probes_agg$adjP = tmp_probes_agg$Score*length(tmp_dm_probes)
tmp_probes_agg$adjP = p.adjust(tmp_probes_agg$adjP, method = "fdr")
aggregated_probes = rownames(filter(tmp_probes_agg, adjP < 0.05))

selectionRes_meth <- list(d_probe_rra = aggregated_probes, gsea_probe_rra = tmp_associated_probes, resampling = T)
saveRDS(selectionRes_meth, file.path(args[[10]], paste0("selectionRes_Methylomics_", args[[9]],".rds")))

