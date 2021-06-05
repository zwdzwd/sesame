#' Add probes to mask
#'
#' This function essentially merge existing probe masking
#' with new prboes to mask
#' 
#' @param sdf a \code{SigDF}
#' @param probes a vector of probe IDs or a logical vector with TRUE
#' representing masked probes
#' @return a \code{SigDF} with added mask
#' @examples
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sum(sdf$mask)
#' sum(addMask(sdf, c("cg14057072", "cg22344912"))$mask)
#' @export
addMask <- function(sdf, probes) {
    if (is.logical(probes)) {
        sdf$mask[probes[sdf$Probe_ID]] = TRUE
    } else {
        sdf$mask[match(probes, sdf$Probe_ID)] = TRUE
    }
    sdf
}

#' Set mask to only the probes specified
#'
#' @param sdf a \code{SigDF}
#' @param probes a vector of probe IDs or a logical vector with TRUE
#' representing masked probes
#' @return a \code{SigDF} with added mask
#' @examples
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sum(sdf$mask)
#' sum(setMask(sdf, "cg14959801")$mask)
#' sum(setMask(sdf, c("cg14057072", "cg22344912"))$mask)
#' @export
setMask <- function(sdf, probes) {
    addMask(resetMask(sdf), probes)
}

#' Reset Masking
#'
#' @param sdf a \code{SigDF}
#' @return a new \code{SigDF} with mask reset to all FALSE
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sum(sdf$mask)
#' sdf <- addMask(sdf, c("cg14057072", "cg22344912"))
#' sum(sdf$mask)
#' sum(resetMask(sdf)$mask)
#' @export
resetMask <- function(sdf) {
    sdf$mask = FALSE
    sdf
}

#' Mask beta values by design quality
#' 
#' Currently quality masking only supports three platforms
#' 
#' @param sdf a \code{SigDF} object
#' @param mask.use.manifest use manifest to mask probes
#' @param mask.use.tcga whether to use TCGA masking, only applies to HM450
#' @return a filtered \code{SigDF}
#' @examples
#' sesameDataCache("EPIC") # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sum(sdf$mask)
#' sum(qualityMask(sdf)$mask)
#' @export 
qualityMask <- function(
    sdf, mask.use.manifest = TRUE,
    manifest = NULL,
    mask.use.tcga = FALSE) {

    ## mask using manifest
    if (mask.use.manifest) {
        if (is.null(manifest)) {
            manifest <- sesameDataGet(paste0(
                platform(sdf), '.address'))$ordering
        }
        if ("mask" %in% colnames(manifest)) {
            sdf <- addMask(sdf, setNames(manifest$mask, manifest$Probe_ID))
        }
    }

    ## mask HM450/HM27/EPIC using TCGA masking
    if (mask.use.tcga) {
        if(!(platform(sdf) %in% c('HM27','HM450','EPIC'))) {
            message(sprintf(
                "TCGA masking is not supported for %s.", platform(sdf)))
            return(sdf)
        }
        stopifnot(platform(sdf) == 'HM450')
        sdf = addMask(sdf, sesameDataGet('HM450.probeInfo')$mask.tcga)
    }

    sdf
}

