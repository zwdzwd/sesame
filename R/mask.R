
#' Reset Masking
#'
#' @param sset a \code{SigSet}
#' @return a new \code{SigSet} with mask reset to empty
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sset.no.mask <- resetMask(sset)
#' @export
resetMask <- function(sset) {
    probes <- probeNames(sset)
    sset@extra$mask <- setNames(rep(FALSE, length(probes)), probes)
    sset
}

#' Save current mask
#'
#' @param sset a \code{SigSet} object
#' @param to new mask name
#' @return a new \code{SigSet} object
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sset <- resetMask(sset)
#' sset <- saveMask(sset)
#' @export
saveMask <- function(sset, to='mask2') {
    ## TODO: change the examples to more sset with meaningful mask
    sset@extra[[to]] <- sset@extra$mask
    sset
}


#' Save current mask
#'
#' @param sset a \code{SigSet} object
#' @param from name of a previously saved mask
#' @return a new \code{SigSet} object
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sset <- resetMask(sset)
#' sset <- saveMask(sset)
#' sset <- restoreMask(sset)
#' @export
restoreMask <- function(sset, from='mask2') {
    ## TODO: change the examples to more sset with meaningful mask
    sset@extra$mask <- sset@extra[[from]]
    sset
}

#' Mask beta values by design quality
#' 
#' Currently quality masking only supports three platforms
#' 
#' @param sset a \code{SigSet} object
#' @param mask.use.tcga whether to use TCGA masking, only applies to HM450
#' @return a filtered \code{SigSet}
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sset.masked <- qualityMask(sset)
#' @export 
qualityMask <- function(
    sset,
    mask.use.tcga = FALSE) {

    if(!(sset@platform %in% c('HM27','HM450','EPIC'))) {
        message(sprintf(
            "Quality masking is not supported for %s.", sset@platform))
        return(sset)
    }
        
    
    if (mask.use.tcga) {
        stopifnot(sset@platform == 'HM450')
        masked <- sesameDataGet('HM450.probeInfo')$mask.tcga
    } else {
        masked <- sesameDataGet(paste0(sset@platform, '.probeInfo'))$mask
    }

    if (!extraHas(sset, 'mask')) {
        resetMask(sset);
    }
    sset@extra$mask[masked] <- TRUE

    return(sset)
}

#' Mask Sigset by detection p-value
#'
#' @param sset a \code{SigSet}
#' @param pval.method which method to use in calculating p-values
#' @param pval.threshold the p-value threshold
#' @return a filtered \code{SigSet}
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sset.masked <- detectionMask(sset)
#' @export
detectionMask <- function(
    sset, pval.method=NULL, pval.threshold=0.05) {
    if (is.null(pval.method)) {
        pv <- pval(sset)
    } else {
        stopifnot(
            extraHas(sset, 'pvals') &&
                pval.method %in% names(sset@extra$pvals))
        pv <- sset@extra$pvals[[pval.method]]
    }

    if (!extraHas(sset, 'mask')) {
        resetMask(sset);
    }
    sset@extra$mask[pv[names(sset@extra$mask)] > pval.threshold] <- TRUE

    sset
}
