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
        sdf$mask[probes[sdf$Probe_ID]] <- TRUE
    } else {
        sdf$mask[match(probes, sdf$Probe_ID)] <- TRUE
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
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sum(sdf$mask)
#' sdf <- addMask(sdf, c("cg14057072", "cg22344912"))
#' sum(sdf$mask)
#' sum(resetMask(sdf)$mask)
#' @export
resetMask <- function(sdf) {
    sdf$mask <- FALSE
    sdf
}

#' list existing quality masks for a SigDF
#'
#' @param sdf a SigDF object
#' @return a tibble of masks
#' @examples
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' listAvailableMasks(sdf)
#' @export
listAvailableMasks <- function(sdf) {
    KYCG_getDBs(sprintf(
        "%s.mask", sdfPlatform(sdf)), summary=TRUE)
}

#' Mask beta values by design quality
#' 
#' Currently quality masking only supports three platforms
#' see also listAvailableMasks(sdf)
#' 
#' @param sdf a \code{SigDF} object
#' @param mask_groups mask group name, default to "recommended"
#' @return a filtered \code{SigDF}
#' @examples
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sum(sdf$mask)
#' sum(qualityMask(sdf)$mask)
#'
#' ## list available masks, the mask_name column
#' listAvailableMasks(sdf)
#' 
#' @export 
qualityMask <- function(sdf, mask_groups = "recommended") {
    masks <- do.call(c, KYCG_getDBs(sprintf(
        "%s.mask", sdfPlatform(sdf)), "recommended", silent=TRUE))
    
    addMask(sdf, masks)
}

