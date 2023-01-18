#' caretLogLoss
#'
#' @param data
#' @param lev
#' @param model
#'
#' @description Selfmade LogLoss function for the package
#'
#' @return
#' @export
#'
#' @examples
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
