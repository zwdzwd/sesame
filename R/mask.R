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
#' @param verbose print more messages
#' @return a new \code{SigDF} with mask reset to all FALSE
#' @examples
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sum(sdf$mask)
#' sdf <- addMask(sdf, c("cg14057072", "cg22344912"))
#' sum(sdf$mask)
#' sum(resetMask(sdf)$mask)
#' @export
resetMask <- function(sdf, verbose = FALSE) {
    sdf$mask <- FALSE
    sdf
}

#' list existing quality masks for a SigDF
#'
#' @param platform EPIC, MM285, HM450 etc
#' @param verbose print more messages
#' @return a tibble of masks
#' @examples
#' listAvailableMasks("EPIC")
#' @export
listAvailableMasks <- function(platform, verbose = FALSE) {
    stopifnot(is.character(platform))
    KYCG_getDBs(sprintf(
        "%s.mask", platform), summary=TRUE, silent=!verbose)
}

## list probes masked by probe nonuniqueness
nonuniqMask <- function(platform, verbose = FALSE) {
    stopifnot(is.character(platform))
    df <- listAvailableMasks(platform, verbose = verbose)
    if(is.null(df)) { return(NULL) }
    mask_names <- c("nonunique", "sub35_copy", "multi", "design_issue")
    mask_names <- df$mask_name[df$mask_name %in% mask_names]
    if (length(mask_names) > 0) {
        do.call(c, KYCG_getDBs(sprintf("%s.mask", platform),
            mask_names, silent=!verbose))
    } else {
        NULL
    }
}

#' Mask beta values by design quality
#' 
#' Currently quality masking only supports three platforms
#' see also listAvailableMasks(sdfPlatform(sdf))
#' 
#' @param sdf a \code{SigDF} object
#' @param mask_names mask names, default to "recommended", can be a vector
#' of multiple masks, e.g., c("design_issue", "multi"), NULL to skip
#' @param prefixes mask by probe ID prefixes, e.g., cg
#' @param verbose print more messages
#' @return a filtered \code{SigDF}
#' @examples
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sum(sdf$mask)
#' sum(qualityMask(sdf)$mask)
#' sum(qualityMask(sdf, mask_names = NULL, prefixes = "rs")$mask)
#'
#' ## list available masks, the mask_name column
#' listAvailableMasks(sdfPlatform(sdf))
#' 
#' @export 
qualityMask <- function(sdf, mask_names = "recommended", prefixes = NULL,
    verbose = FALSE) {

    masks <- character(0)
    if (!is.null(prefixes)) { # mask by prefix
        masks <- c(masks, do.call(c, lapply(prefixes, function(pfx) {
            grep(sprintf("^%s", pfx), sdf$Probe_ID, value=TRUE) })))
    }

    ## mask by predefined sets
    platform <- sdfPlatform(sdf, verbose = verbose)
    df <- listAvailableMasks(platform, verbose = verbose)
    if(is.null(df)) { return(sdf) }
    mask_names <- df$mask_name[df$mask_name %in% mask_names]
    if (length(mask_names) > 0 && !is.null(mask_names)) {
        masks <- c(masks, 
            do.call(c, KYCG_getDBs(sprintf(
                "%s.mask", platform),
                mask_names, silent=!verbose)))
    }
    
    addMask(sdf, masks)
}

#' Mask SigDF by probe ID prefix
#'
#' @param sdf SigDF
#' @param prefixes prefix characters
#' @param invert use the complement set
#' @return SigDF
#' @examples
#' sdf <- resetMask(sesameDataGet("MM285.1.SigDF"))
#' sum(prefixMask(sdf, c("ctl","rs"))$mask)
#' sum(prefixMask(sdf, c("ctl"))$mask)
#' sum(prefixMask(sdf, c("ctl","rs","ch"))$mask)
#' @export
prefixMask <- function(sdf, prefixes = NULL, invert = FALSE) {
    idx <- Reduce("|", lapply(prefixes, function(pfx) {
        grepl(sprintf("^%s", pfx), sdf$Probe_ID) }))
    
    if (invert) {
        sdf[!idx, "mask"] <- TRUE
    } else {
        sdf[idx, "mask"] <- TRUE
    }
    sdf
}

#' Mask all but C probes in SigDF
#'
#' @param sdf SigDF
#' @return SigDF
#' @examples
#' sdf <- resetMask(sesameDataGet("MM285.1.SigDF"))
#' sum(prefixMaskButC(sdf)$mask)
#' @export
prefixMaskButC <- function(sdf) {
    prefixMask(sdf, c("cg", "ch"), invert = TRUE)
}

#' Mask all but CG probes in SigDF
#'
#' @param sdf SigDF
#' @return SigDF
#' @examples
#' sdf <- resetMask(sesameDataGet("MM285.1.SigDF"))
#' sum(prefixMaskButCG(sdf)$mask)
#' @export
prefixMaskButCG <- function(sdf) {
    prefixMask(sdf, "cg", invert = TRUE)
}
