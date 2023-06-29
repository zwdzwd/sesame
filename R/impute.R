
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
