#' Convert human arrays to previous platforms
#' Missing probes are replaced using NAs.
#'
#' @param sdf SigDF data frame
#' @param platform HM450 or EPIC
#' @return a new SigDF for the older platform
#' @examples
#' sdf <- sesameDataGet("EPIC.5.SigDF.normal")[[1]]
#' sdf_out <- convertTo(sdf, "HM450")
#' @export
convertTo <- function(sdf, platform=c("HM450", "EPIC")) {
    probes <- sdf$Probe_ID
    if (sdfPlatform(sdf) == "EPICv2") {
        probes <- vapply(strsplit(probes, "_"),
            function(x) x[1], character(1))
    }

    if (platform == "HM450") {
        sdf_ref <- sesameDataGet("HM450.10.SigDF")[[1]]
    } else if (platform == "EPIC") {
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
#' @param platform platform of target imputation
#' @param tissue tissue context of imputation
#' @param probes Probe ID if not set to vector names
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
imputeTo <- function(betas, platform=NULL, probes=NULL, tissue="Blood") {
    if (is.null(probes)) {
        if (is.matrix(betas)) {
            probes <- rownames(betas)
        } else {
            probes <- names(betas)
        }
    }
    platform <- sesameData_check_platform(platform, probes = probes)

    if (platform %in% c("HM450", "EPIC") &&
        any(grepl("_", grep("^cg", probes, value=TRUE)))) {
        betas <- betasCollapseToPfx(betas)
        if (is.matrix(betas)) {
            probes <- rownames(betas)
        } else {
            probes <- names(betas)
        }
    }
    
    ddefault <- sesameDataGet(sprintf("%s.imputationDefault", platform))
    ddefault <- ddefault[[tissue]]
    dat <- setNames(ddefault[["median"]], ddefault$Probe_ID)
    idx <- match(ddefault$Probe_ID, probes)
    dat[!is.na(idx)] <- betas[na.omit(idx)]
    dat
}
