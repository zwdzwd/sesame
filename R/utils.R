#' Convert beta-value to M-value
#'
#' Convert beta-value to M-value (aka logit transform)
#' @param b vector of beta values
#' @return a vector of M values
#' @export
BetaValueToMValue <- function(b) {
  log2(b/(1-b))
}

#' Convert M-value to beta-value
#'
#' Convert M-value to beta-value (aka inverse logit transform)
#' @param m a vector of M values
#' @return a vector of beta values
#' @export
MValueToBetaValue <- function(m) {
  2^m/(1+2^m)
}

## print message
smessage <- function(...) {
  cat('[', as.character(Sys.time()),'] ', ..., '\n', sep='')
}
