#' a List of SigSets with some methods of its own
#' 
#' @import S4Vectors
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
  l <- List(...) 
  platform <- unique(sapply(l, slot, "platform"))
  if (length(platform) > 1) {
    stop("All SigSet elements in a SigSetList must be from the same platform")
  }
  new("SigSetList", l, platform=platform, elementType="SigSet")
}


#' read an entire directory's worth of IDATs into a SigSetList 
#' 
#' @param   path      the path from which to read IDATs (default ".")
#' @param   parallel  run in parallel? (default FALSE) 
#' 
#' @return            a SigSetList 
#'
#' @export
SigSetListFromPath <- function(path=".", parallel=FALSE) {
  idats <- list.files(path=path, patt="idat")
  stubs <- unique(sub(".gz", "", sub("_(Grn|Red).idat", "", idats)))
  names(stubs) <- stubs
  if (parallel == TRUE) {
    SigSetList(mclapply(stubs, readIDATpair, verbose=TRUE))
  } else { 
    SigSetList(lapply(stubs, readIDATpair, verbose=TRUE))
  }
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
            cat("platform: ", slot(object, "platform"), "\n")
            # left open for later use
          })
