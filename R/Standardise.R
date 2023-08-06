#' Standardise
#'
#' @description Standardises a data frame
#'
#' @param dataFrame to be standardised
#'
#' @return a standardised data frame
#'
#' @author Ulrich Asemann
Standardise <- function(dataFrame){
  # Remove columns where sd() is 0
  remove <- apply(as.data.frame(dataFrame), 2, function(x) sd(x) == 0)
  dataFrame <- dataFrame[,!remove]

  # Standardise the data frame values
  dataFrame <- apply(as.data.frame(dataFrame), 2, function(x) (x - mean(x))/sd(x)) %>% as.matrix()
  return(dataFrame)
}
