
#' Impute to platform
#' 
#' @param betas named vector of beta values
#' @param platform platform of target imputation
#' @param tissue tissue context of imputation
#' @param probes Probe ID if not set to vector names
#' @examples
#'
#' betas <- c("cg04707299"=0.2, "cg13380562"=0.9, "cg00000103"=0.1)
#' betas2 <- imputeTo(betas, "EPIC")
#' betas2[names(betas)]
#' head(betas2)
imputeTo <- function(betas, platform=NULL, probes=NULL, tissue="Blood") {
    if (is.null(probes)) {
        probes <- names(betas)
    }
    platform <- sesameData_check_platform(platform, probes = probes)
    ddefault <- sesameDataGet(sprintf("%s.imputationDefault", platform))
    ddefault <- ddefault[[tissue]]
    dat <- setNames(ddefault[["median"]], ddefault$Probe_ID)
    idx <- match(ddefault$Probe_ID, probes)
    dat[!is.na(idx)] <- betas[na.omit(idx)]
    dat
}
