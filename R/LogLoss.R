#' LogLoss Function
#'
#' @description The self written LogLoss function for the package.
#'
#' @param pred
#' @param true
#' @param eps
#' @param weights
#'
#' @return Returns the LogLoss.
#'
#' @author Ulrich Asemann

LogLoss <- function(pred,
                    true,
                    eps = 1e-15,
                    weights = NULL) {
  pred = pmin(pmax(pred, eps), 1 - eps) # Bound the results
  if (is.null(weights)) {
    return(-(sum(
      true * log(pred) + (1 - true) * log(1 - pred)
    )) / length(true))
  } else{
    return(-weighted.mean(true * log(pred) + (1 - true) * log(1 - pred), weights))
  }
}
