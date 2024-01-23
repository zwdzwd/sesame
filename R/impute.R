
#' lift over beta values or SigDFs to another Infinium platform
#' This function wraps ID conversion and provide optional
#' imputation functionality.
#'
#' @param x either named beta value (vector or matrix), probe IDs
#' or SigDF(s)
#' if input is a matrix, probe IDs should be in the row names
#' if input is a numeric vector, probe IDs should be in the vector
#' names.
#' If input is a character vector, the input will be
#' considered probe IDs.
#' @param target_platform the platform to take the data to
#' @param source_platform optional information of the source data
#' platform (when there might be ambiguity).
#' @param mapping a liftOver mapping file. Typically this file
#' contains empirical evidence whether a probe mapping is reliable.
#' If given, probe ID-based mapping will be skipped. This is to
#' perform more stringent probe ID mapping.
#' @param impute whether to impute or not, default is FALSE
#' @param celltype the cell type / tissue context of imputation,
#' if not given, will use nearest neighbor to find out.
#' @return imputed data, vector, matrix, SigDF(s)
#' @examples
#'
#' \dontrun{
#' sesameDataCache()
#' sdf = sesameDataGet("EPICv2.8.SigDF")[["GM12878_206909630042_R08C01"]]
#' betas = openSesame(sdf)
#' betas_HM450 = liftOver(betas, "HM450", impute=TRUE)
#'
#' ## directly map probes
#' cg_epic2 = grep("cg", names(sesameData_getManifestGRanges("EPICv2")), value=T)
#' head(liftOver(cg_epic2, "HM450"))
#'
#' cg_hm450 = grep("cg", names(sesameData_getManifestGRanges("HM450")), value=T)
#' head(liftOver(cg_hm450, "EPICv2"))
#'
#' rs_epic2 = grep("rs", names(sesameData_getManifestGRanges("EPICv2")), value=T)
#' head(liftOver(rs_epic2, "HM450", source_platform="EPICv2"))
#' 
#' }
#' @export
liftOver <- function(x,
    target_platform, source_platform=NULL,
    mapping=NULL, impute=FALSE, celltype="Blood") {
    
    if (is.numeric(x) || is.matrix(x)) {
        imputeTo(x, target_platform, mapping=mapping, impute=impute,
            celltype=celltype)
    } else if (is(x, "SigDF")) {
        convertTo(x, target_platform)
    } else if (is.character(x)) {
        source_platform <- sesameData_check_platform(source_platform, x)
        y <- names(sesameData_getManifestGRanges(target_platform))
        if (target_platform %in% c("EPIC", "HM450", "HM27") &&
            source_platform %in% c("EPICv2", "MSA")) {
            x2 <- vapply(strsplit(x, "_"), function(xx) xx[1], character(1))
            y2 <- y
        } else if (target_platform %in% c("EPICv2", "MSA") &&
                   source_platform %in% c("EPIC", "HM450", "HM27")) {
            x2 <- x
            y2 <- vapply(strsplit(y, "_"), function(xx) xx[1], character(1))
        } else {
            x2 <- x
            y2 <- y
        }
        idx <- match(y2, x2)
        y <- y[!is.na(idx)]
        names(y) <- x[idx[!is.na(idx)]]
        y
    } else if (is.list(x) && is(x[[1]], "SigDF")) {
        lapply(x, convertTo, platform = target_platform)
    }
}

#' Convert human arrays to previous platforms
#' Missing probes are replaced using NAs.
#'
#' @param sdf SigDF data frame
#' @param target_platform HM450 or EPIC
#' @return a new SigDF for the older platform
#' @examples
#' sdf <- sesameDataGet("EPIC.5.SigDF.normal")[[1]]
#' sdf_out <- convertTo(sdf, "HM450")
#' @export
convertTo <- function(sdf, target_platform=c("HM450", "EPIC")) {
    probes <- sdf$Probe_ID
    if (sdfPlatform(sdf) == "EPICv2") {
        probes <- vapply(strsplit(probes, "_"),
            function(x) x[1], character(1))
    }

    if (target_platform == "HM450") {
        sdf_ref <- sesameDataGet("HM450.10.SigDF")[[1]]
    } else if (target_platform == "EPIC") {
        sdf_ref <- sesameDataGet("EPIC.5.SigDF.normal")[[1]]
    }

    idx <- match(sdf_ref$Probe_ID, probes)
    sdf_out <- sdf_ref
    sdf_out$MG <- sdf$MG[idx]
    sdf_out$MR <- sdf$MR[idx]
    sdf_out$UG <- sdf$UG[idx]
    sdf_out$UR <- sdf$UR[idx]
    sdf_out$col <- sdf$col[idx]
    sdf_out$mask <- sdf$mask[idx]
    sdf_out$mask[is.na(sdf_out$mask)] <- TRUE
    sdf_out
}

#' Impute to platform
#' 
#' @param betas named vector or matrix of beta values
#' @param target_platform platform of target imputation
#' @param mapping a liftOver mapping file, if given, probe ID-based
#' mapping will be skipped.
#' @param celltype celltype/tissue context of imputation, if not given, will
#' use nearest neighbor to determine.
#' @param impute whether to impute or not, default is no
#' @return imputed data, vector or matrix
#' @examples
#'
#' betas <- c("cg04707299"=0.2, "cg13380562"=0.9, "cg00000103"=0.1)
#' betas_imputed <- imputeTo(betas, "HM450")
#' 
#' betas <- setNames(seq(0,1,length.out=3),
#'     c("cg00004963_TC21", "cg00004963_TC22", "cg00004747_TC21"))
#' betas_imputed <- imputeTo(betas, "HM450")
#' 
#' @export
imputeTo <- function(betas, target_platform=NULL, mapping=NULL,
    impute=FALSE, celltype="Blood") {

    btm <- betas
    if (!is.matrix(btm)) { btm <- cbind(btm) }
    target_platform <- sesameData_check_platform(target_platform)

    if (!is.null(mapping)) { # use liftOver file
        btm <- btm[mapping$ID_source,,drop=FALSE]
        rownames(btm) <- mapping$ID_target
    } else if ( # convert by mapping prefixes
        target_platform %in% c("HM450", "EPIC") &&
        any(grepl("_", grep("^cg", rownames(btm), value=TRUE)))) {
        btm <- betasCollapseToPfx(btm)
    }

    probes2 <- sesameDataGet(sprintf(
        "%s.address", target_platform))$ordering$Probe_ID
    btm <- btm[match(probes2, rownames(btm)),,drop=FALSE]
    rownames(btm) <- probes2
    
    if (impute) {
        dd0 <- sesameDataGet(sprintf(
            "%s.imputationDefault", target_platform))
        dd <- dd0$data[[celltype]]
        idx <- match(dd0$Probe_ID, rownames(btm))
        btm <- apply(btm[idx,,drop=FALSE], 2, function(x) {
            idx1 <- (is.na(idx) | is.na(x))
            x[idx1] <- dd$median[idx1]
            x })
        rownames(btm) <- dd0$Probe_ID
    }

    if (!is.matrix(betas)) { btm <- btm[,1] }
    btm
}
