#' Add probes to mask
#'
#' This function essentially merge existing probe masking
#' with new probes to mask
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
    KYCG_getDBs(sprintf("%s.Mask", platform),
        summary=TRUE, silent=!verbose, ignore.case = TRUE)$dbname
}

## list probes masked by probe nonuniqueness
nonuniqMask <- function(platform, verbose = FALSE) {
    stopifnot(is.character(platform))
    dbnames <- listAvailableMasks(platform, verbose = verbose)
    if(is.null(dbnames)) { return(NULL) }
    mask_names <- c("M_nonuniq", "nonunique",
        "sub35_copy", "multi", "design_issue")
    mask_names <- dbnames[dbnames %in% mask_names]
    if (length(mask_names) > 0) {
        do.call(c, KYCG_getDBs(sprintf("%s.Mask", platform),
            mask_names, silent=!verbose))
    } else {
        NULL
    }
}

#' Recommended mask names for each Infinium platform
#'
#' The returned name is the db name used in KYCG.mask
#' @return a named list of mask names
#' @examples
#' recommendedMaskNames()[["EPIC"]]
#' recommendedMaskNames()[["EPICv2"]]
#' 
#' @export
recommendedMaskNames <- function() {
    list(
        EPICv2 = c(
            "M_1baseSwitchSNPcommon_5pt",
            "M_2extBase_SNPcommon_5pt",
            "M_mapping", "M_nonuniq", "M_SNPcommon_5pt"),
        MM285 = c("ref_issue", "nonunique", "design_issue"),
        EPIC = c(
            "mapping", "channel_switch", "snp5_GMAF1p",
            "extension", "sub30_copy"),
        HM450 = c(
            "mapping", "channel_switch", "snp5_GMAF1p",
            "extension", "sub30_copy"),
        HM27 = c("mask"))
}

#' get probe masking by mask names
#' 
#' @param platform EPICv2, EPIC, HM450, HM27, ...
#' @param mask_names mask names (see listAvailableMasks), by default: NULL
#' Note that setting this does not turn off recommended masking
#' to turn off recommend masking, you need
#' mask_names = "<your mask names>", use_recommended = FALSE
#' to turn off all masking, you need
#' mask_names = NULL, use_recommended = FALSE
#' @param use_recommended whether or not to apply recommmended masking,
#' by default: TRUE.
#' see recommendedMaskNames() for detail.
#' @return a vector of probe ID
#' @examples
#'
#' recommendedMaskNames()[["EPIC"]]
#' length(getMask("EPIC"))
#' length(getMask("HM450"))
#' length(getMask("MM285"))
#' 
#' @export
getMask <- function(platform = "EPICv2",
    mask_names = NULL, use_recommended = TRUE) {
    stopifnot(is.character(platform))
    masks <- NULL
    if (use_recommended) {
        masks <- c(masks, KYCG_getDBs(sprintf("%s.Mask", platform),
            recommendedMaskNames()[[platform]], silent=TRUE, ignore.case=TRUE))
    }

    if (length(mask_names) > 0) {
        masks <- c(masks, KYCG_getDBs(sprintf("%s.Mask", platform),
            mask_names, silent=TRUE, ignore.case=TRUE))
    }

    if (is.null(masks)) {
        NULL
    } else {
        unique(do.call(c, masks))
    }
}

#' Mask beta values by design quality
#' 
#' Currently quality masking only supports three platforms
#' see also listAvailableMasks(sdfPlatform(sdf))
#' 
#' @param sdf a \code{SigDF} object
#' @param verbose print more messages
#' @param ... masking details see getMask()
#' @return a filtered \code{SigDF}
#' @examples
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' sum(sdf$mask)
#' sum(qualityMask(sdf)$mask)
#' sum(qualityMask(sdf, mask_names = NULL)$mask)
#'
#' ## list available masks, the dbname column
#' listAvailableMasks(sdfPlatform(sdf))
#' listAvailableMasks("EPICv2")
#' 
#' @export 
qualityMask <- function(sdf, verbose = FALSE, ...) {

    ## mask by predefined sets
    platform <- sdfPlatform(sdf, verbose = verbose)
    masks <- getMask(platform, ...)
    if (is.null(masks)) {
        return(sdf)
    } else {
        addMask(sdf, masks)
    }
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
