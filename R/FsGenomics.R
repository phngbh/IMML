#' Feature Selection Genomics Data
#'
#' @description Feature selection for genomics data using GSEA.
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
#' @param phenotype The name of the column in phenotypeIDs, which will be used
#'   in the analysis.
#' @param binaryFilePrefix The file path to the .bim/.bed files used. Only
#'   specify the prefix without the file types.
#' @param geneLocFile The file path to the gene location file to be used.
#' @param geneSetFile The file path to the gene set file to be used.
#' @param outputPrefix The prefix of the output file.
#' @param seed The seed used for random number generation. Using the same seed
#'   ensures reproducibility.
#'
#' @return A text file will be saved in the environment holding the results.
#'
#' @author Wilhelm Glaas
#'
#' @export

FsGenomics <- function(trainIDs,
                       testIDs,
                       dataIDs,
                       phenotypeIDs,
                       phenotype,
                       binaryFilePrefix,
                       geneLocFile,
                       geneSetFile,
                       outputPrefix,
                       seed = 123) {

  # make samples dir
  dir.create('samples')
  for (i in 1:length(trainIDs)) {
    write.table(data.frame(V1 = trainIDs[[i]], V2 = trainIDs[[i]]),
                file = file.path(getwd(), paste0("samples/samples_",i,".txt")),
                col.names = F, row.names = F, sep = "\t")
  }
  allIds = c(trainIDs[[1]], testIDs[[1]])
  write.table(data.frame(V1 = allIds, V2 = allIds),
              file = file.path(getwd(), "samples/samples_fs.txt"),
              col.names = F, row.names = F, sep = "\t")

  # make phenotype file
  merged = merge(rownames_to_column(dataIDs, var = 'ID'),
                 rownames_to_column(phenotypeIDs, var = 'ID'), by = "ID")
  merged = data.frame(merged$Genomics, merged$Genomics, merged[phenotype]) %>%
            drop_na()
  colnames(merged) = c('FID', 'IID', 'CKD')
  merged$CKD = merged$CKD + 1
  merged = merged[order(merged[['CKD']], merged[['FID']]),]
  write.table(merged, 'phenotype_file.txt', quote = F, row.names = F, sep = '\t')
  rm(merged)

  # run magma script
  command = paste(
    system.file('magma.sh', package = 'IMLpackage'),
    getwd(),
    binaryFilePrefix,
    geneLocFile,
    file.path(getwd(), 'samples'),
    file.path(getwd(), 'phenotype_file.txt'),
    geneSetFile,
    outputPrefix,
    system.file(package = 'IMLpackage'),
    length(trainIDs)
  )
  system(command)

  return("Feature selection genomics done!")
}
