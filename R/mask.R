#' Add probes to mask
#'
#' This function essentially merge existing probe masking
#' with new prboes to mask
#' 
#' @param sset a \code{SigSet}
#' @param probes a vector of probe IDs or a logical vector with TRUE
#' representing masked probes
#' @return a \code{SigSet} with added mask
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sum(sesame::mask(sset))
#' sum(sesame::mask(addMask(sset, c("cg14057072", "cg22344912"))))
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

## all new versions should run this while initializing sset
initializeMask <- function(sset) {
    if (is.null(sset@extra[["mask"]])) {
        pnames <- probeNames(sset)
        sset@extra$mask <- setNames(rep(FALSE, length(pnames)), pnames)
    }
    sset
}

#' Set mask to only the probes specified
#'
#' @param sset a \code{SigSet}
#' @param probes a vector of probe IDs or a logical vector with TRUE
#' representing masked probes
#' @return a \code{SigSet} with added mask
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
#' @return a new \code{SigSet} with mask reset to all FALSE
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sum(mask(sset))
#' sset <- addMask(sset, c("cg14057072", "cg22344912"))
#' sum(mask(sset))
#' sum(mask(resetMask(sset)))
#' @export
resetMask <- function(sset) {
    pnames <- probeNames(sset)
    sset@extra$mask <- setNames(rep(FALSE, length(pnames)), pnames)
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
#' @param manifest the manifest file that contains mask column
#' @param mask.use.tcga whether to use TCGA masking, only applies to HM450
#' @return a filtered \code{SigSet}
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' sum(mask(sset))
#' sset.masked <- qualityMask(sset)
#' sum(mask(sset.masked))
#' @export 
qualityMask <- function(
    sset,
    mask.use.manifest = TRUE,
    manifest = NULL,
    mask.use.tcga = FALSE) {

    ## shouldn't have to run this for lastest sigsets
    if (!extraHas(sset, 'mask') || length(sset@extra$mask) == 0) {
        sset <- resetMask(sset);
    }

    ## mask using manifest
    if (mask.use.manifest) {
        if (is.null(manifest)) {
            manifest <- sesameDataGet(paste0(
                sset@platform, '.address'))$ordering
        }
        if ("mask" %in% colnames(manifest)) {
            sset <- addMask(sset, setNames(manifest$mask, manifest$Probe_ID))
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

