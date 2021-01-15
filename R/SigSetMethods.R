#' @include sesame.R

#' probeNames generic
#'
#' @param x object of \code{Sigset}
#'
#' @rdname probeNames-methods
#' @docType methods
#' @export
setGeneric("probeNames", function(x) {
    standardGeneric("probeNames")
})

#' Get Probe Names of SigSet class
#'
#' @rdname probeNames-methods
#' @return A char vector
#' @aliases probeNames,SigSet-method
#' @examples
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' head(probeNames(sset))
setMethod("probeNames", "SigSet", function(x) {
    c(rownames(x@II), rownames(x@IR), rownames(x@IG))
})

##############
#### pval ####
##############

#' pval replacement generic
#' 
#' @param x object of \code{SigSet}
#' @param value new value
#' 
#' @rdname pval-replace-methods
#' @docType methods
#' @export
setGeneric("pval<-", function(x, value) {
    standardGeneric("pval<-")
})

#' Replace pval slot of SigSet class
#' 
#' @rdname pval-replace-methods
#' @return a new \code{SigSet}
#' @docType methods
#' @aliases pval<-,SigSet-method
#' @examples 
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' df <- pval(sset)
#' df[1] <- 0.01
#' pval(sset) <- df
setReplaceMethod("pval", "SigSet", function(x, value){ x@pval <- value; x})

#' pval getter generic
#' 
#' @param x object of \code{SigSet}
#' 
#' @rdname pval-methods
#' @docType methods
#' @export
setGeneric("pval", function(x) {
    standardGeneric("pval")
})

#' Get pval slot of SigSet class
#' 
#' @rdname pval-methods
#' @return The pval slot of \code{SigSet}
#' @aliases pval,SigSet-method
#' @examples 
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' head(pval(sset))
setMethod("pval", "SigSet", function(x){ x@pval })

############
#### IG ####
############
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
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' df <- IG(sset)
#' df[1,1] <- 10
#' IG(sset) <- df 
setReplaceMethod("IG", "SigSet", function(x, value){ x@IG <- value; x})

#' IG getter generic
#' 
#' @param x object of \code{SigSet}
#' 
#' @rdname IG-methods
#' @docType methods
#' @export
setGeneric("IG", function(x) {
    standardGeneric("IG")
})

#' Get IG slot of SigSet class
#' 
#' @rdname IG-methods
#' @return The IG slot of \code{SigSet}
#' @aliases IG,SigSet-method
#' @examples 
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' head(IG(sset))
setMethod("IG", "SigSet", function(x){ x@IG })

#' Get IG slot of SigSet class that passes mask
#' This is the same as IG() if there is no mask set
#'
#' @param sset SigSet object
#' @return The IG slot that passes extra$mask filter
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' head(IGpass(sset))
#' @export
IGpass <- function(sset) {
    stopifnot(is(sset, "SigSet"))
    if (extraHas(sset, "mask")) {
        sset@IG[!extraGet(sset,"mask")[rownames(sset@IG)],]
    } else {
        sset@IG
    }
}

############
#### IR ####
############
#' IR replacement generic
#' 
#' @param x object of \code{SigSet}
#' @param value new value
#' 
#' @rdname IR-replace-methods
#' @docType methods
#' @export
setGeneric("IR<-", function(x, value) {
    standardGeneric("IR<-")
})

#' Replace IR slot of SigSet class
#' 
#' @rdname IR-replace-methods
#' @return a new \code{SigSet}
#' @docType methods
#' @aliases IR<-,SigSet-method
#' @examples 
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' df <- IR(sset)
#' df[1,1] <- 10
#' IR(sset) <- df 
setReplaceMethod("IR", "SigSet", function(x, value){ x@IR <- value; x})

#' IR getter generic
#' 
#' @param x object of \code{SigSet}
#' 
#' @rdname IR-methods
#' @docType methods
#' @export
setGeneric("IR", function(x) {
    standardGeneric("IR")
})

#' Get IR slot of SigSet class
#' 
#' @rdname IR-methods
#' @return The IR slot of \code{SigSet}
#' @aliases IR,SigSet-method
#' @examples 
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' head(IR(sset))
setMethod("IR", "SigSet", function(x){ x@IR })

#' Get IR slot of SigSet class that passes mask
#' This is the same as IR() if there is no mask set
#'
#' @param sset SigSet object
#' @return The IR slot that passes extra$mask filter
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' head(IRpass(sset))
#' @export
IRpass <- function(sset) {
    stopifnot(is(sset, "SigSet"))
    if (extraHas(sset, "mask")) {
        sset@IR[!extraGet(sset,"mask")[rownames(sset@IR)],]
    } else {
        sset@IR
    }
}

############
#### II ####
############
#' II replacement generic
#' 
#' @param x object of \code{SigSet}
#' @param value new value
#' 
#' @rdname II-replace-methods
#' @docType methods
#' @export
setGeneric("II<-", function(x, value) {
    standardGeneric("II<-")
})

#' Replace II slot of SigSet class
#' 
#' @rdname II-replace-methods
#' @return a new \code{SigSet}
#' @docType methods
#' @aliases II<-,SigSet-method
#' @examples 
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' df <- II(sset)
#' df[1,1] <- 10
#' II(sset) <- df 
setReplaceMethod("II", "SigSet", function(x, value){ x@II <- value; x})

#' II getter generic
#' 
#' @param x object of \code{SigSet}
#' 
#' @rdname II-methods
#' @docType methods
#' @export
setGeneric("II", function(x) {
    standardGeneric("II")
})

#' Get II slot of SigSet class
#' 
#' @rdname II-methods
#' @return The II slot of \code{SigSet}
#' @aliases II,SigSet-method
#' @examples 
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' head(II(sset))
setMethod("II", "SigSet", function(x){ x@II })

#' Get II slot of SigSet class that passes mask
#' This is the same as II() if there is no mask set
#'
#' @param sset SigSet object
#' @return The II slot that passes extra$mask filter
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' head(IIpass(sset))
#' @export
IIpass <- function(sset) {
    stopifnot(is(sset, "SigSet"))
    if (extraHas(sset, "mask")) {
        sset@II[!extraGet(sset,"mask")[rownames(sset@II)],]
    } else {
        sset@II
    }
}

##############
#### oobG ####
##############
#' oobG replacement generic
#' 
#' @param x object of \code{SigSet}
#' @param value new value
#' 
#' @rdname oobG-replace-methods
#' @docType methods
#' @export
setGeneric("oobG<-", function(x, value) {
    standardGeneric("oobG<-")
})

#' Replace oobG slot of SigSet class
#' 
#' @rdname oobG-replace-methods
#' @return a new \code{SigSet}
#' @docType methods
#' @aliases oobG<-,SigSet-method
#' @examples 
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' df <- oobG(sset)
#' df[1,1] <- 10
#' oobG(sset) <- df 
setReplaceMethod("oobG", "SigSet", function(x, value){ x@oobG <- value; x})

#' oobG getter generic
#' 
#' @param x object of \code{SigSet}
#' 
#' @rdname oobG-methods
#' @docType methods
#' @export
setGeneric("oobG", function(x) {
    standardGeneric("oobG")
})

#' Get oobG slot of SigSet class
#' 
#' @rdname oobG-methods
#' @return The oobG slot of \code{SigSet}
#' @aliases oobG,SigSet-method
#' @examples 
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' head(oobG(sset))
setMethod("oobG", "SigSet", function(x){ x@oobG })

#' Get oobG slot of SigSet class that passes mask
#' This is the same as oobG() if there is no mask set
#'
#' @param sset SigSet object
#' @return The oobG slot that passes extra$mask filter
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' head(oobGpass(sset))
#' @export
oobGpass <- function(sset) {
    stopifnot(is(sset, "SigSet"))
    if (extraHas(sset, "mask")) {
        sset@oobG[!extraGet(sset,"mask")[rownames(sset@oobG)],]
    } else {
        sset@oobG
    }
}

##############
#### oobR ####
##############
#' oobR replacement generic
#' 
#' @param x object of \code{SigSet}
#' @param value new value
#' 
#' @rdname oobR-replace-methods
#' @docType methods
#' @export
setGeneric("oobR<-", function(x, value) {
    standardGeneric("oobR<-")
})

#' Replace oobR slot of SigSet class
#' 
#' @rdname oobR-replace-methods
#' @return a new \code{SigSet}
#' @docType methods
#' @aliases oobR<-,SigSet-method
#' @examples 
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' df <- oobR(sset)
#' df[1,1] <- 10
#' oobR(sset) <- df 
setReplaceMethod("oobR", "SigSet", function(x, value){ x@oobR <- value; x})

#' oobR getter generic
#' 
#' @param x object of \code{SigSet}
#' 
#' @rdname oobR-methods
#' @docType methods
#' @export
setGeneric("oobR", function(x) {
    standardGeneric("oobR")
})

#' Get oobR slot of SigSet class
#' 
#' @rdname oobR-methods
#' @return The oobR slot of \code{SigSet}
#' @aliases oobR,SigSet-method
#' @examples 
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' head(oobR(sset))
setMethod("oobR", "SigSet", function(x){ x@oobR })

#' Get oobR slot of SigSet class that passes mask
#' This is the same as oobR() if there is no mask set
#'
#' @param sset SigSet object
#' @return The oobR slot that passes extra$mask filter
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' head(oobRpass(sset))
#' @export
oobRpass <- function(sset) {
    stopifnot(is(sset, "SigSet"))
    if (extraHas(sset, "mask")) {
        sset@oobR[!extraGet(sset,"mask")[rownames(sset@oobR)],]
    } else {
        sset@oobR
    }
}

############
#### ctl ###
############
#' ctl replacement generic
#' 
#' @param x object of \code{SigSet}
#' @param value new value
#' 
#' @rdname ctl-replace-methods
#' @docType methods
#' @export
setGeneric("ctl<-", function(x, value) {
    standardGeneric("ctl<-")
})

#' Replace ctl slot of SigSet class
#' 
#' @rdname ctl-replace-methods
#' @return a new \code{SigSet}
#' @docType methods
#' @aliases ctl<-,SigSet-method
#' @examples 
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' df <- ctl(sset)
#' df[1,1] <- 10
#' ctl(sset) <- df 
setReplaceMethod("ctl", "SigSet", function(x, value){ x@ctl <- value; x})

#' ctl getter generic
#' 
#' @param x object of \code{SigSet}
#' 
#' @rdname ctl-methods
#' @docType methods
#' @export
setGeneric("ctl", function(x) {
    standardGeneric("ctl")
})

#' Get ctl slot of SigSet class
#' 
#' @rdname ctl-methods
#' @return The ctl slot of \code{SigSet}
#' @aliases ctl,SigSet-method
#' @examples 
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' head(ctl(sset))
setMethod("ctl", "SigSet", function(x){ x@ctl })

##############
#### extra ###
##############
#' extra replacement generic
#' 
#' @param x object of \code{SigSet}
#' @param value new value
#' 
#' @rdname extra-replace-methods
#' @docType methods
#' @export
setGeneric("extra<-", function(x, value) {
    standardGeneric("extra<-")
})

#' Replace extra slot of SigSet class
#' 
#' @rdname extra-replace-methods
#' @return a new \code{SigSet}
#' @docType methods
#' @aliases extra<-,SigSet-method
#' @examples 
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' df <- extra(sset)
#' extra(sset) <- list(pval=numeric(0))
setReplaceMethod("extra", "SigSet", function(x, value){ x@extra <- value; x})

#' extra getter generic
#' 
#' @param x object of \code{SigSet}
#' 
#' @rdname extra-methods
#' @docType methods
#' @export
setGeneric("extra", function(x) {
    standardGeneric("extra")
})

#' Get extra slot of SigSet class
#' 
#' @rdname extra-methods
#' @return The extra slot of \code{SigSet}
#' @aliases extra,SigSet-method
#' @examples 
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' head(extra(sset))
setMethod("extra", "SigSet", function(x){ x@extra })

