
#' lift over beta values or SigDFs to another Infinium platform
#' This function wraps ID conversion and provide optional
#' imputation functionality.
#' 
#' @param x either beta value (vector or matrix) or SigDF(s)
#' @param target_platform the platform to take the data to
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
#' sesameDataCache()
#' ## sdf = sesameDataGet("EPICv2.8.SigDF")[["GM12878_206909630042_R08C01"]]
#' ## betas = openSesame(sdf)
#' ## betas_HM450 = liftOver(betas, "HM450", impute=TRUE)
#' @export
liftOver <- function(x,
    target_platform, mapping=NULL, impute=FALSE, celltype="Blood") {
    
    if (is.numeric(x) || is.matrix(x)) {
        imputeTo(x, target_platform, mapping=mapping, impute=impute,
            celltype=celltype)
    } else if (is(x, "SigDF")) {
        convertTo(x, target_platform)
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
    
    if (!is.matrix(betas)) {
        betas <- cbind(betas)
    }
    target_platform <- sesameData_check_platform(target_platform)

    if (!is.null(mapping)) { # use liftOver file
        betas <- betas[mapping$ID_source,,drop=FALSE]
        rownames(betas) <- mapping$ID_target
    } else if ( # convert by mapping prefixes
        target_platform %in% c("HM450", "EPIC") &&
        any(grepl("_", grep("^cg", rownames(betas), value=TRUE)))) {
        betas <- betasCollapseToPfx(betas)
    }

    probes2 <- sesameDataGet(sprintf(
        "%s.address", target_platform))$ordering$Probe_ID
    betas <- betas[match(probes2, rownames(betas)),,drop=FALSE]
    rownames(betas) <- probes2
    
    if (impute) {
        dd0 <- sesameDataGet(sprintf(
            "%s.imputationDefault", target_platform))
        dd <- dd0$data[[celltype]]
        idx <- match(dd0$Probe_ID, rownames(betas))
        betas <- apply(betas[idx,,drop=FALSE], 2, function(x) {
            idx1 <- (is.na(idx) | is.na(x))
            x[idx1] <- dd$median[idx1]
            x })
        rownames(betas) <- dd0$Probe_ID
    }
    betas
}
