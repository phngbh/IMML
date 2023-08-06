#' BedToDataframeCopyNumber
#'
#' @description
#'
#' @param bed filepath to the ".bed" file
#'
#' @return
#'
#' @author Ulrich Asemann
BedToDataframeCopyNumber <- function(bed){
  # Load data from path
  data <- read.plink(bed = bed)
  dataFrame <- as.data.frame(data$genotypes)
  # Get rownames
  rownames <- rownames(dataFrame)
  dataFrame.char <- apply(dataFrame,2,as.numeric)
  rownames(dataFrame.char) <- rownames
  # Change values
  for (i in 1:nrow(dataFrame.char)){
    for (j in 1:ncol(dataFrame.char)){
      if (dataFrame.char[i,j] == 1){
        dataFrame.char[i,j] = 2
      } else if (dataFrame.char[i,j] == 2){
        dataFrame.char[i,j] = 1
      } else if (dataFrame.char[i,j] == 3){
        dataFrame.char[i,j] = 0
      }
    }
  }
  return(dataFrame.char)
}
