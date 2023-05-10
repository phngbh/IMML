#' caretLogLoss Function
#'
#' @description The selfwritten caretLogLoss function for the package.
#'
#' @param data
#' @param lev
#' @param model
#'
#' @return Returns the caretLogLoss.
#'
#' @author Ulrich Asemann

caretLogLoss <- function(data, lev = NULL, model = NULL) {
  cls <- levels(data$obs) #find class names
  loss <- LogLoss(
    pred = data[, cls[2]],
    true = as.numeric(data$obs) - 1,
    weights = data$weights
  )
  names(loss) <- c('myLogLoss')
  loss
}
