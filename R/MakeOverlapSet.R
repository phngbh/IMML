#' MakeOverlapSet
#'
#' @description
#'
#' @param dataList
#' @param IDTable
#' @param phenotypeIDs
#'
#' @return
#'
#' @author Ulrich Asemann

MakeOverlapSet <- function(dataList,
                           IDTable,
                           phenotypeIDs) {
  listFull <-
    na.omit(IDTable) %>% filter(as.character(Clinical) %in% rownames(phenotypeIDs)) %>% lapply(., as.character)
  IDFil <- list()
  for (i in names(listFull)) {
    IDFil[[i]] <- which(listFull[[i]] %in% rownames(dataList[[i]]))
  }
  # IDFil = lapply(listFull, function(x, dataList = dataList) which(x %in% rownames(dataList[[x]])))
  IDFinal <- Reduce(intersect, IDFil)
  for (i in names(dataList)) {
    dataList[[i]] <- dataList[[i]][listFull[[i]][IDFinal],]
    rownames(dataList[[i]]) <-
      rownames(dataList[["Clinical"]][listFull[["Clinical"]][IDFinal],])
  }
  # dataList <-
  #   lapply(dataList, function(x, IDFinal = IDFinal)
  #     return(x[IDFinal, ]))
  return(dataList)
}
