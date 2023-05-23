# code to use feature selection methods
fs_transcriptomics(
  training_IDs,
  testing_IDs,
  sample_IDs,
  phenotypes,
  transcriptomics_processed,
  trans_gene_information,
  trans_lowest_level_pathways
)

fs_metabolomics(
  training_IDs,
  testing_IDs,
  sample_IDs,
  phenotypes,
  metabolomics_processed,
  meta_pathway_list,
  meta_lowest_level_pathways,
  meta_anno_fil
)

fs_clinical(training_IDs,
            testing_IDs,
            sample_IDs,
            phenotypes,
            clinical_processed)

fs_proteomics(training_IDs,
              testing_IDs,
              sample_IDs,
              phenotypes,
              proteomics_processed,
              prot_lowest_level_pathways,
              prot_somamer_info_edited)

fs_methylomics(training_IDs,
               testing_IDs,
               sample_IDs,
               phenotypes,
               methylomics_processed,
               meth_lowest_level_pathways)

fs_genomics(geno_lowest_level_pathways,
            geno_annot)

# clear enviroment
rm(list = ls())

# save enviroment
save.image(file = 'metabolomicsEnv.RData')

# load csv file from raw github
library(readr)
urlfile <-
  "https://raw.githubusercontent.com/phngbh/DSPN/master/3.FeatureSelection/Proteomics/somamer_info_edited.csv"
mydata <- read.csv(url(urlfile), head = TRUE, sep = ";")
write.csv(mydata, "somamer_info_edited.csv", row.names = FALSE)
read.csv("somamer_info_edited.csv")

# Testing of functions and other stuff
# Will be deleted later on


resamples = training_IDs$`Feature Selection IDs`$Transcriptomics
samples = unlist(resamples) %>% unique()
phenotypes[as.character(sample_IDs$Clinical)[match(samples, as.character(sample_IDs$Transcriptomics))], "inc3", drop = F]

x = (c(1, 2, 2, 2, 3, 4, 4))
y = (c(2, 3, 4, 4, 4, 5))
intersect(x, y)

sd <- 0.3 * sqrt(4 / rchisq(100, df = 4))
y <- matrix(rnorm(100 * 6, sd = sd), 100, 6)
rownames(y) <- paste("Gene", 1:100)
y[1:2, 4:6] <- y[1:2, 4:6] + 2
design <- cbind(Grp1 = 1, Grp2vs1 = c(0, 0, 0, 1, 1, 1))
options(digits = 3)
lmFit(y, design)

drop_na(as.data.frame(phenotypes$inc3))


samples <-
  training_IDs$`Feature Selection IDs`$Transcriptomics %>% unlist() %>% unique()

info <-
  phenotypes[as.character(sample_IDs$Clinical)[match(samples, as.character(sample_IDs$Transcriptomics))],
             "inc3", drop = F]

info

test <- kable(gsea_final)

# test fs_metabolomics
toppw <-
  fs_metabolomics(
    training_IDs,
    testing_IDs,
    sample_IDs,
    phenotypes,
    metabolomics_processed,
    meta_pathway_list,
    meta_lowest_level_pathways,
    meta_anno_fil
  )
topde <-
  fs_metabolomics(
    training_IDs,
    testing_IDs,
    sample_IDs,
    phenotypes,
    metabolomics_processed,
    meta_pathway_list,
    meta_lowest_level_pathways,
    meta_anno_fil
  )

ranklist <- topde$t
topde$ChEBI <-
  meta_anno_fil$ChEBI[match(topde$ID, meta_anno_fil$Mnumber)]
names(ranklist) <- topde$ChEBI
ranklist <- sort(ranklist)

geneset_reactome <- reactomePathways(names(ranklist))
geneset_reactome <-
  geneset_reactome[intersect(names(geneset_reactome), toppw)]

fgseaRes <- fgsea(
  pathways = geneset_reactome,
  stats    = ranklist,
  minSize  = 1,
  maxSize  = 200
) %>% arrange(pval)

fgseaRes <-
  fgseaRes %>% mutate(
    leadingEdge = mapIdsList(
      x = org.Hs.eg.db,
      keys = leadingEdge,
      keytype = "ENTREZID",
      column = "SYMBOL"
    )
  )

edge = fgseaRes$leadingEdge %>% unlist() %>% unique() %>% na.omit()
probe = meta_anno_fil$Mnumber[meta_anno_fil$ChEBI %in% edge]

# clinical feature selection


# testing parallel computing

# DESCRIPTION
foreach,
parallel,
doSNOW

library(foreach)
library(doSNOW)

detectCores()


cl <- makeCluster(2)
registerDoSNOW(cl)

foreach(i=1:1000) %dopar% {
  i
}

stopCluster(cl)


x = c("1","2")









