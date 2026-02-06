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

KYCG_listDBGroups <- function(filter = NULL, path = NULL, type = NULL) {

    if (is.null(path)) {
        gps <- sesameDataList("KYCG", full=TRUE)[,c("Title","Description")]
        gps$type <- vapply(strsplit(
            gps$Description, " "), function(x) x[2], character(1))
        gps$Description <- str_replace(
            gps$Description, "KYCG categorical database holding ", "")
        if (!is.null(filter)) {
            gps <- gps[grepl(filter, gps$Title),]
        }
        if (!is.null(type)) {
            gps <- gps[gps$type %in% type,]
        }
    } else {
        gps <- basename(list.files(path, recursive = TRUE))
    }
    gps
}

guess_dbnames <- function(nms, platform = NULL,
    allow_multi = FALSE, type = NULL, silent = FALSE,
    ignore.case = FALSE) {
    
    gps <- KYCG_listDBGroups(type = type)
    nms <- do.call(c, lapply(nms, function(nm) {
        if (nm %in% gps$Title) {
            return(nm)
        } else if (length(grep(nm, gps$Title, ignore.case=ignore.case)) >= 1) {
            ret <- grep(nm, gps$Title, value=TRUE, ignore.case=ignore.case)
            if (!allow_multi) { ret <- ret[1]; }
            return(ret)
        } else if (length(grep(nm, gps$Title, ignore.case=ignore.case)) == 0) {
            res <- gps$Title[apply(do.call(cbind, lapply(
                strsplit(nm, "\\.")[[1]],
                function(q1) grepl(q1, gps$Title, ignore.case=ignore.case))),
                1, all)]
            if (length(res) == 1) {
                return(res[1])
            }
        }
        return(nm)
    }))
    if (!is.null(platform)) {
        nms <- grep(platform, nms, value = TRUE)
    }
    if (!silent) {
        message("Selected the following database groups:")
        invisible(lapply(seq_along(nms), function(i) {
            message(sprintf("%d. %s", i, nms[i]))
        }))
    }
    nms
}

## now an internal function
KYCG_getDBs <- function(group_nms, db_names = NULL, platform = NULL,
    summary = FALSE, allow_multi = FALSE,
    ignore.case = FALSE, type = NULL, silent = FALSE) {
    
    if (!is.character(group_nms)) {
        return(group_nms)
    }

    group_nms <- guess_dbnames(group_nms, platform = platform,
        allow_multi = TRUE, type = type, silent = silent,
        ignore.case = ignore.case)
    ## group_nms <- group_nms[sesameDataHas(group_nms)]
    group_nms <- group_nms[group_nms %in% sesameDataList()$Title]
    if (length(group_nms) == 0) {
        return(NULL)
    }
    res <- do.call(c, lapply(unname(group_nms), function(nm) {
        dbs <- sesameDataGet(nm)
        setNames(lapply(seq_along(dbs), function(ii) {
            db <- dbs[[ii]]
            attr(db, "group") <- nm
            attr(db, "dbname") <- names(dbs)[ii]
            db
        }), names(dbs))}))
    
    if (summary) {
        do.call(bind_rows, lapply(res, attributes))
    } else if (is.null(db_names)) {
        res
    } else {
        stopifnot(all(db_names %in% names(res)))
        res[db_names]
    }
}

#' list existing quality masks for a SigDF
#'
#' @param platform EPIC, MM285, HM450 etc
#' @param verbose print more messages
#' @return a tibble of masks
#' @examples
#' listAvailableMasks("EPICv2")
#' @importFrom reshape2 melt
#' @export
listAvailableMasks <- function(platform, verbose = FALSE) {
    stopifnot(is.character(platform))
    KYCG_getDBs(sprintf("%s.Mask", platform),
        summary=TRUE, silent=!verbose, ignore.case = TRUE)$dbname
}

## these are probes excluded from background subtraction and pOOBAH etc.
backgroundMask <- function(platform, verbose = FALSE) {
    stopifnot(is.character(platform))
    dbnames <- listAvailableMasks(platform, verbose = verbose)
    if(is.null(dbnames)) { return(NULL) }
    mask_names <- c("M_nonuniq", "nonunique",
        "sub35_copy", "multi", "design_issue",
        "M_1baseSwitchSNPcommon_5pt",
        "M_1baseSwitchSNPcommon_1pt")
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
#' @param platform string for platform, e.g., EPIC, EPICv2, MSA
#' @return a named list of mask names
#' @examples
#' recommendedMaskNames("EPICv2")
#' recommendedMaskNames("EPIC")
#' 
#' @export
recommendedMaskNames <- function(platform) {
    ma <- list(
        MM285 = c("ref_issue", "nonunique", "design_issue"),
        EPIC = c(
            "mapping", "channel_switch", "snp5_GMAF1p",
            "extension", "sub30_copy"),
        EPICv2 = c(
            "M_1baseSwitchSNPcommon_5pt",
            "M_2extBase_SNPcommon_5pt",
            "M_mapping", "M_nonuniq", "M_SNPcommon_5pt"),
        MSA = c(
            "M_1baseSwitchSNPcommon_5pt",
            "M_2extBase_SNPcommon_5pt",
            "M_mapping", "M_nonuniq", "M_SNPcommon_5pt"),
        HM450 = c(
            "mapping", "channel_switch", "snp5_GMAF1p",
            "extension", "sub30_copy"),
        HM27 = c("mask"))
    if (platform %in% names(ma)) {
        ma[[platform]]
    } else {
        ma[["EPICv2"]]
    }
}

#' get probe masking by mask names
#' 
#' @param platform EPICv2, EPIC, HM450, HM27, ...
#' @param mask_names mask names (see listAvailableMasks)
#' by default: "recommended"
#' see recommendedMaskNames() for detail.
#' @return a vector of probe ID
#' @examples
#'
#' length(getMask("MSA", "recommended"))
#' length(getMask("EPICv2", "recommended"))
#' length(getMask("EPICv2", c("recommended", "M_SNPcommon_1pt")))
#' length(getMask("EPICv2", "M_mapping"))
#' length(getMask("EPIC"))
#' length(getMask("HM450"))
#' length(getMask("MM285"))
#' 
#' @export
getMask <- function(platform = "EPICv2", mask_names = "recommended") {
    stopifnot(is.character(platform))
    res <- lapply(mask_names, function(mask_name) {
        if (mask_name == "recommended") {
            res <- KYCG_getDBs(sprintf("%s.Mask", platform),
                recommendedMaskNames(platform),
                silent=TRUE, ignore.case=TRUE)
        } else {
            res <- KYCG_getDBs(sprintf("%s.Mask", platform),
                mask_name, silent=TRUE, ignore.case=TRUE)
        }
        if (!is.null(res)) {
            res <- do.call(c, res)
        }
        res
    })
    if (!is.null(res)) {
        res <- unique(do.call(c, res))
    }
    res
}

#' Mask beta values by design quality
#' 
#' Currently quality masking only supports three platforms
#' see also listAvailableMasks(sdfPlatform(sdf))
#' 
#' @param sdf a \code{SigDF} object
#' @param mask_names a vector of masking groups, see listAvailableMasks
#' use "recommended" for recommended masking. One can also combine
#' "recommended" with other masking groups by specifying a vector, e.g.,
#' c("recommended", "M_mapping")
#' @param verbose be verbose
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
qualityMask <- function(sdf,
    mask_names="recommended", verbose=TRUE) {
    
    ## mask by predefined sets
    platform <- sdfPlatform(sdf, verbose=verbose)
    masks <- getMask(platform, mask_names=mask_names)
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
