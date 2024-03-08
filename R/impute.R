#' Impute of missing data of specific platform
#' 
#' @param betas named vector of beta values
#' @param platform platform
#' @param celltype celltype/tissue context of imputation, if not given, will
#' use nearest neighbor to determine.
#' @param sd_max maximum standard deviation in imputation confidence
#' @return imputed data, vector or matrix
imputeBetas <- function(betas, platform = NULL,
    celltype = NULL, sd_max = 999) {

    df <- sesameDataGet(sprintf("%s.imputationDefault", platform))
    d2q <- match(names(betas), df$Probe_ID)
    celltype <- names(which.max(vapply(df$data, function(x) cor(
        betas, x$median[d2q], use="na.or.complete"), numeric(1))))
    if (is.null(celltype)) {
        celltype <- "Blood"
    }
    idx <- is.na(betas)
    mn <- df$data[[celltype]]$median[d2q][idx]
    sd <- df$data[[celltype]]$sd[d2q][idx]
    mn[sd > sd_max] <- NA
    betas[idx] <- mn
    betas
}
