#' Add probes to mask. This function essentially merge existing probe masking
#' with new prboes to mask
#' 
#' @param sset a \code{SigSet}
#' @param probes a vector of probe IDs or a logical vector with TRUE
#' representing masked probes
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sum(mask(sset))
#' sum(mask(addMask(sset, c("cg14057072", "cg22344912"))))
#' @export
addMask <- function(sset, probes) {
    if (!extraHas(sset, "mask") || length(sset@extra$mask) == 0) {
        sset <- resetMask(sset)
    }
    if (is.logical(probes)) {
        sset@extra$mask[probes[match(probeNames(sset), names(probes))]] = TRUE
    } else {
        sset@extra$mask[match(probes, probeNames(sset))] = TRUE
    }
    sset
}

#' Set mask to only the probes specified
#'
#' @param sset a \code{SigSet}
#' @param probes a vector of probe IDs or a logical vector with TRUE
#' representing masked probes
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sum(mask(sset))
#' sum(mask(setMask(sset, "cg14959801")))
#' sum(mask(setMask(sset, c("cg14057072", "cg22344912"))))
#' @export
setMask <- function(sset, probes) {
    addMask(resetMask(sset), probes)
}

#' Reset Masking
#'
#' @param sset a \code{SigSet}
#' @return a new \code{SigSet} with mask reset to empty
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sum(mask(sset))
#' sset <- addMask(sset, c("cg14057072", "cg22344912"))
#' sum(mask(sset))
#' sum(mask(resetMask(sset)))
#' @export
resetMask <- function(sset) {
    sset@extra$mask <- rep(FALSE, length(probeNames(sset)))
    sset
}

#' Save current mask
#'
#' @param sset a \code{SigSet} object
#' @param to new mask name
#' @return a new \code{SigSet} object
#' @examples
#' sesameDataCache("EPIC") # if not done yet
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
#' sesameDataCache("EPIC") # if not done yet
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
#' @param mask.use.manifest use manifest to mask probes
#' @param mask.use.tcga whether to use TCGA masking, only applies to HM450
#' @return a filtered \code{SigSet}
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sset.masked <- qualityMask(sset)
#' @export 
qualityMask <- function(
    sset,
    mask.use.manifest = TRUE,
    mask.use.tcga = FALSE) {

    if (!extraHas(sset, 'mask') || length(sset@extra$mask) == 0) {
        sset <- resetMask(sset);
    }

    ## mask using manifest
    if (mask.use.manifest && extraHas(sset, "maskManifest")) {
        if (!extraHas(sset, "mask") || length(sset@extra$mask) == 0) {
            sset@extra$mask = sset@extra$maskManifest
        } else {
            mask = sset@extra$maskManifest[names(sset@extra$mask)]
            sset@extra$mask[mask] = TRUE
        }
    }

    ## mask HM450/HM27/EPIC using TCGA masking
    if (mask.use.tcga) {
        if(!(sset@platform %in% c('HM27','HM450','EPIC'))) {
            message(sprintf(
                "TCGA masking is not supported for %s.", sset@platform))
            return(sset)
        }
        stopifnot(sset@platform == 'HM450')
        masked <- sesameDataGet('HM450.probeInfo')$mask.tcga
        sset@extra$mask[masked] = TRUE
    }

    return(sset)
}

