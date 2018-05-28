#' IG replacement generic
#' 
#' @param x object of \code{SigSet}
#' @param value new value
#' 
#' @rdname IG-replace-methods
#' @docType methods
#' @export
setGeneric("IG<-", function(x, value) {
    standardGeneric("IG<-")
})

#' Replace IG slot of SigSet class
#' 
#' @rdname IG-replace-methods
#' @return a new \code{SigSet}
#' @docType methods
#' @aliases IG<-,SigSet-method
#' @examples 
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')
#' df <- IG(sset)
#' df[1,1] <- 10
#' IG(sset) <- df 
setReplaceMethod("IG", "SigSet", function(x, value){ x@IG <- value})

#' IG getter generic
#' 
#' @param x object of \code{SigSet}
#' 
#' @rdname IG-methods
#' @docType methods
#' @export
setGeneric("IG", function(x) {
    standardGeneric("IG")
    x@IG
})

#' Get IG slot of SigSet class
#' 
#' @rdname IG-methods
#' @return The IG slot of \code{SigSet}
#' @aliases IG,SigSet-method
#' @examples 
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')
#' head(IG(sset))
setMethod("IG", "SigSet", function(x){ x@IG })

