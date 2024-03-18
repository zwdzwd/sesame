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

#' Impute Missing Values with Mean
#' This function replaces missing values (NA) in a matrix, default is row
#' means.
#' 
#' @param mx A matrix
#' @param axis A single integer. Use 1 to impute column means (default),
#' and 2 to impute row means.
#' @return A matrix with missing values imputed.
#' @examples
#' mx <- cbind(c(1, 2, NA, 4), c(NA, 2, 3, 4))
#' imputeBetasMatrixByMean(mx, axis = 1)
#' imputeBetasMatrixByMean(mx, axis = 2)
#' @export
imputeBetasMatrixByMean <- function(mx, axis = 1) {
    stopifnot(is.matrix(mx))
    if (axis == 1) {
        t(apply(mx, 1, function(x) {
            x[is.na(x)] <- mean(x, na.rm = TRUE);
            x
        }))
    } else if (axis == 2) {
        apply(mx, 2, function(x) {
            x[is.na(x)] <- mean(x, na.rm = TRUE);
            x
        })
    } else {
        stop("Invalid axis. Use 1 for columns or 2 for rows.")
    }
}
