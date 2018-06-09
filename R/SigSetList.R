#' a List of SigSets with some methods of its own
#' 
#' @import S4Vectors
#' @import Matrix
#' 
#' @exportClass SigSetList
setClass("SigSetList",
         representation(platform="character"),
         contains="SimpleList")


#' constructor
#' 
#' @param   ...       the SigSet objects that will be the List elements
#' 
#' @return            a SigSetList 
#'
#' @export
SigSetList <- function(...) { 
  platform <- unique(sapply(list(...), slot, "platform"))
  if (length(platform) > 1) {
    stop("All SigSet elements in a SigSetList must be from the same platform")
  }
  new("SigSetList", 
      List(...),
      platform=platform,
      elementType="SigSet")
}


#' SigSetList methods (centralized).
#'  
#' `show`         print a summary of the SigSetList.
#'
#' @param x       a SigSetList 
#' @param sset    a SigSetList 
#' @param object  a SigSetList 
#'
#' @name  SigSetList-methods
NULL


#' @rdname SigSetList-methods
#' @export
setMethod("show", signature(object="SigSetList"),
          function(object) {
            callNextMethod()
            masked <- apply(slot(object, "mask"), 2, is.na)
            maskedpct <- paste0(round(masked/nrow(slot(object,"mask"))), "%")
            cat(S4Vectors:::labeledLine("masked", masked))
          })
