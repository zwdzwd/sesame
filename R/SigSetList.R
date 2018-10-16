#' a List of SigSets with some methods of its own
#' 
#' @importFrom S4Vectors List
#' 
#' @exportClass SigSetList
setClass(
    "SigSetList",
    representation(platform="character"),
    contains="SimpleList")


#' constructor
#' 
#' @param ... the SigSet objects that will be the List elements
#' @return a SigSetList 
#' @examples
#' sset1 <- readIDATpair(file.path(system.file(
#'     'extdata','',package='sesameData'), '4207113116_A'))
#' 
#' sset2 <- readIDATpair(file.path(system.file(
#'     'extdata','',package='sesameData'), '4207113116_B'))
#' 
#' SigSetList(sset1, sset2)
#' @export
SigSetList <- function(...) {
    l <- List(...)
    platform <- unique(vapply(l, slot, character(1), "platform"))
    if (length(platform) > 1) {
        stop("All SigSetList elements must be from the same platform")
    }
    new("SigSetList", l, platform=platform, elementType="SigSet")
}


#' read an entire directory's worth of IDATs into a SigSetList 
#' 
#' @param   path      the path from which to read IDATs (default ".")
#' @param   recursive whether to search recursively
#' @param   parallel  run in parallel? (default FALSE) 
#' 
#' @return            a SigSetList
#' @examples
#' ## Load all IDATs from directory
#' ssets <- SigSetListFromPath(
#'     system.file("extdata", "", package = "sesameData"))
#' @export
SigSetListFromPath <- function(path=".", parallel=FALSE, recursive=TRUE) {
    # idats <- list.files(path=path, pattern="idat")
    # stubs <- unique(sub(".gz", "", sub("_(Grn|Red).idat", "", idats)))
    # names(stubs) <- stubs
    stubs <- searchIDATprefixes(path, recursive = recursive)
    message("Found ", length(stubs), " IDAT files in ", path, ".")
    SigSetListFromIDATs(stubs=stubs, parallel=parallel)
}


#' read IDATs into a SigSetList 
#' 
#' FIXME: switch from `parallel` to BiocParallel
#' 
#' @param   stubs       the IDAT filename stubs
#' @param   parallel    run in parallel? (default FALSE) 
#' 
#' @return              a SigSetList 
#'
#' @importFrom parallel mclapply
#' @examples 
#' ## a SigSetList of length 1
#' ssets <- SigSetListFromIDATs(file.path(
#'     system.file("extdata", "", package = "sesameData"), "4207113116_A"))
#' @export
SigSetListFromIDATs <- function(stubs, parallel=FALSE) {
    if (parallel == TRUE) {
        SigSetList(mclapply(stubs, readIDATpair, verbose=TRUE))
    } else { 
        SigSetList(lapply(stubs, readIDATpair, verbose=TRUE))
    }
}


#' SigSetList methods (centralized). 
#' Currently scarce...
#'  
#' `show`         print a summary of the SigSetList.
#'
#' @param object  a SigSetList 
#'
#' @name SigSetList-methods
NULL

#' @rdname SigSetList-methods
#' @return Description of SigSetList
#' @examples 
#' SigSetListFromPath(system.file("extdata", "", package = "sesameData"))
#' @export
setMethod(
    "show", signature(object="SigSetList"),
    function(object) {
        callNextMethod()
        platform <- slot(object, "platform")
        probes <- switch(
            platform,
            "EPIC"=865918,
            "HM450"=485577,
            "HM27"=27578)
        cat("platform:", platform, paste0("(", probes, " probes)"), "\n")
    })
