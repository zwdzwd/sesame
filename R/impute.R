#' Impute of missing data of specific platform
#' 
#' @param betas named vector of beta values
#' @param platform platform
#' @param celltype celltype/tissue context of imputation, if not given, will
#' use nearest neighbor to determine.
#' @param sd_max maximum standard deviation in imputation confidence
#' @param BPPARAM use MulticoreParam(n) for parallel processing
#' @return imputed data, vector or matrix
#' @examples
#' betas = openSesame(sesameDataGet("EPIC.1.SigDF"))
#' sum(is.na(betas))
#' betas2 = imputeBetas(betas, "EPIC")
#' sum(is.na(betas2))
#' 
#' @export
imputeBetas <- function(betas, platform = NULL, BPPARAM = SerialParam(),
    celltype = NULL, sd_max = 999) {

    if (is.matrix(betas)) {
        cnms <- colnames(betas)
        betas <- do.call(cbind, bplapply(seq_len(ncol(betas)), function(i) {
            imputeBetas(betas[,i], platform = platform,
                celltype = celltype, sd_max = sd_max)}, BPPARAM=BPPARAM))
        colnames(betas) <- cnms
        return(betas)
    }

    platform <- sesameData_check_platform(platform, names(betas))
    df <- sesameDataGet(sprintf("%s.imputationDefault", platform))
    d2q <- match(names(betas), df$Probe_ID)
    if (is.null(celltype)) {
        celltype <- names(which.max(vapply(df$data, function(x) cor(
            betas, x$median[d2q], use="na.or.complete"), numeric(1))))
    }
    if (length(celltype) == 0 || is.null(celltype)) {
        celltype <- "Blood"
    }
    idx <- is.na(betas)
    mn <- df$data[[celltype]]$median[d2q][idx]
    sd <- df$data[[celltype]]$sd[d2q][idx]
    mn[sd > sd_max] <- NA
    betas[idx] <- mn
    betas
}

#' Impute missing data based on genomic neighbors.
#'
#' @param betas named vector of beta values
#' @param platform platform
#' @param max_neighbors maximum neighbors to use for dense regions
#' @param max_dist maximum distance to count as neighbor
#' @param BPPARAM use MulticoreParam(n) for parallel processing
#' @return imputed data, vector or matrix
#' @importFrom GenomicRanges resize
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @importFrom dplyr summarize
#' @examples
#' betas = openSesame(sesameDataGet("EPICv2.8.SigDF")[[1]])
#' sum(is.na(betas))
#' betas2 = imputeBetasByGenomicNeighbors(betas, "EPICv2")
#' sum(is.na(betas2))
#' 
#' @export
imputeBetasByGenomicNeighbors <- function(betas, platform = NULL,
    BPPARAM = SerialParam(), max_neighbors = 3, max_dist = 10000) {

    platform <- sesameData_check_platform(platform, names(betas))
    mft <- sesameData_getManifestGRanges(platform)
    mft_missing <- mft[names(mft) %in% names(which(is.na(betas)))]
    mft_nonmiss <- mft[names(which(!is.na(betas)))]
    index <- findOverlaps(resize(mft_missing, max_dist), mft_nonmiss)
    gm <- mft_missing[queryHits(index)]
    gn <- mft_nonmiss[subjectHits(index)]
    df <- tibble(
        cg = names(gm), beg_m = start(gm), end_m = end(gm),
        cg_n = names(gn), beg_n = start(gn), end_n = end(gn))

    df$d1 <- df$beg_m - df$end_n - 1
    df$d2 <- df$beg_n - df$end_m - 1
    df$betas <- betas[df$cg_n]
    df$dist <- pmax(df$d1, df$d2)
    df <- summarize(slice_min(group_by(df, .data[['cg']]),
        n = max_neighbors, order_by = .data[['dist']]),
        mbetas = mean(.data[['betas']]))
    betas[df$cg] <- df$mbetas
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
