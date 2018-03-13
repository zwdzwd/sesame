#' Resolve name to object
#' 
#' @param x the target
#' @param ... extra options
#' @return resolved grid object
Resolve <- function(x, ...) {
    UseMethod('Resolve')
}

#' Calculate Text Bounding
#'
#' Calculate bounding box including texts.
#'
#' W.R.T lower left corner of the view port in the unit of points
#' @param x object
#' @param ... extra options
#' @return bounding box coordinates
CalcTextBounding <- function(x, ...) {
    UseMethod('CalcTextBounding', x)
}
