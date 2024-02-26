#' Lift over beta values or SigDFs to another Infinium platform
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
#' @param ... extra arguments, see ?convertProbeID
#' @return imputed data, vector, matrix, SigDF(s)
#' @examples
#'
#' \dontrun{
#' sesameDataCache()
#' 
#' ## lift SigDF
#' 
#' sdf = sesameDataGet("EPICv2.8.SigDF")[["GM12878_206909630042_R08C01"]]
#' dim(liftOver(sdf, "EPICv2"))
#' dim(liftOver(sdf, "EPIC"))
#' dim(liftOver(sdf, "HM450"))
#'
#' sdfs = sesameDataGet("EPICv2.8.SigDF")[1:2]
#' sdfs_hm450 = liftOver(sdfs, "HM450")
#'
#' sdf = sesameDataGet("EPIC.5.SigDF.normal")[[1]]
#' dim(liftOver(sdf, "EPICv2"))
#' dim(liftOver(sdf, "EPIC"))
#' dim(liftOver(sdf, "HM450"))
#'
#' sdf = sesameDataGet("HM450.10.SigDF")[[1]]
#' dim(liftOver(sdf, "EPICv2"))
#' dim(liftOver(sdf, "EPIC"))
#' dim(liftOver(sdf, "HM450"))
#' 
#' ## lift beta values
#'
#' betas = openSesame(sesameDataGet("EPICv2.8.SigDF")[[1]])
#' betas_hm450 = liftOver(betas, "HM450", impute=TRUE)
#' length(betas_hm450)
#' sum(is.na(betas_hm450))
#' betas_hm450 <- liftOver(betas, "HM450", impute=FALSE)
#' length(betas_hm450)
#' sum(is.na(betas_hm450))
#' betas_epic1 <- liftOver(betas, "EPIC", impute=TRUE)
#' length(betas_epic1)
#' sum(is.na(betas_epic1))
#' betas_epic1 <- liftOver(betas, "EPIC", impute=FALSE)
#' length(betas_epic1)
#' sum(is.na(betas_epic1))
#'
#' betas_matrix = openSesame(sesameDataGet("EPICv2.8.SigDF")[1:4])
#' dim(betas_matrix)
#' betas_matrix_hm450 = liftOver(betas_matrix, "HM450", impute=T)
#' dim(betas_matrix_hm450)
#'
#' ## use empirical evidence in liftOver
#' mapping = sesameDataGet("liftOver.EPICv2ToEPIC")
#' betas_matrix = openSesame(sesameDataGet("EPICv2.8.SigDF")[1:4])
#' dim(liftOver(betas_matrix, "EPIC", mapping = mapping))
#' ## compare to without using empirical evidence
#' dim(liftOver(betas_matrix, "EPIC"))
#' 
#' betas <- c("cg04707299"=0.2, "cg13380562"=0.9, "cg00000103"=0.1)
#' head(liftOver(betas, "HM450", impute=TRUE))
#' 
#' betas <- c("cg00004963_TC21"=0, "cg00004963_TC22"=0.5, "cg00004747_TC21"=1.0)
#' betas_hm450 <- liftOver(betas, "HM450", impute=TRUE)
#' head(na.omit(liftOver(betas, "HM450", impute=FALSE)))
#'
#' ## lift probe IDs
#'
#' cg_epic2 = names(sesameData_getManifestGRanges("EPICv2"))
#' head(liftOver(cg_epic2, "HM450"))
#' 
#' cg_epic2 = grep("cg", names(sesameData_getManifestGRanges("EPICv2")), value=T)
#' head(liftOver(cg_epic2, "HM450"))
#'
#' cg_hm450 = grep("cg", names(sesameData_getManifestGRanges("HM450")), value=T)
#' head(liftOver(cg_hm450, "EPICv2"))
#'
#' rs_epic2 = grep("rs", names(sesameData_getManifestGRanges("EPICv2")), value=T)
#' head(liftOver(rs_epic2, "HM450", source_platform="EPICv2"))
#'
#' probes_epic2 = names(sesameData_getManifestGRanges("EPICv2"))
#' head(liftOver(probes_epic2, "EPIC"))
#' head(liftOver(probes_epic2, "EPIC", target_uniq = TRUE))
#' head(liftOver(probes_epic2, "EPIC", include_new = FALSE))
#' head(liftOver(probes_epic2, "EPIC", include_old = FALSE))
#' head(liftOver(probes_epic2, "EPIC", return_mapping=TRUE))
#' 
#' }
#' @export
liftOver <- function(x,
    target_platform, source_platform=NULL,
    mapping=NULL, impute=FALSE, celltype="Blood", ...) {

    if (is.numeric(x)) {
        if (is.matrix(x)) {
            betas <- apply(x, 2, function(xx) liftOver(
                xx, target_platform, source_platform = NULL,
                mapping = mapping, impute = impute, celltype = celltype))
        } else {
            mapping <- convertProbeID(names(x), target_platform,
                mapping = mapping, return_mapping = TRUE, include_new = TRUE)
            betas <- setNames(x[mapping$ID_source], mapping$ID_target)
            if (impute) {
                betas <- imputeBetas(betas, target_platform, celltype = celltype)
            }
        }
        betas
    } else if (is(x, "SigDF")) {
        mapping <- convertProbeID(
            x$Probe_ID, target_platform, return_mapping = TRUE,
            target_uniq = TRUE, include_new = TRUE)
        x2 <- x[match(mapping$ID_source, x$Probe_ID),]
        x2$Probe_ID <- mapping$ID_target
        x2 <- x2[order(x2$Probe_ID),]
        x2$mask[is.na(x2$mask)] <- TRUE
        rownames(x2) <- NULL
        x2
    } else if (is.character(x)) {
        convertProbeID(x, target_platform, source_platform, ...)
    } else if (is.list(x) && is(x[[1]], "SigDF")) {
        lapply(x, liftOver, target_platform, source_platform = source_platform,
            mapping = mapping, impute = impute, celltype = celltype, ...)
    }
}

#' Convert Probe ID
#' 
#' @param x source probe IDs
#' @param target_platform the platform to take the data to
#' @param source_platform optional source platform
#' @param mapping a liftOver mapping file. Typically this file
#' contains empirical evidence whether a probe mapping is reliable.
#' If given, probe ID-based mapping will be skipped. This is to
#' perform more stringent probe ID mapping.
#' @param target_uniq whether the target Probe ID should be kept unique.
#' @param include_new if true, include mapping of added probes
#' @param include_old if true, include mapping of deleted probes
#' @param return_mapping return mapping table, instead of the target IDs.
#' @return mapped probe IDs, or mapping table if return_mapping = T
#' @importFrom dplyr full_join
#' @importFrom dplyr distinct
#' @importFrom stats setNames
convertProbeID <- function(
    x, target_platform, source_platform = NULL, mapping = NULL,
    target_uniq = TRUE, include_new = FALSE, include_old = FALSE,
    return_mapping = FALSE) {

    if (is.null(mapping)) {
        source_platform <- sesameData_check_platform(source_platform, x)
        dfs <- tibble(ID_source = x)
        dft <- tibble(
            ID_target = sesameDataGet(sprintf(
                "%s.address", target_platform))$ordering$Probe_ID)
        if (target_platform %in% c("EPIC", "HM450", "HM27") &&
            source_platform %in% c("EPICv2", "MSA")) {
            dfs$prefix <- vapply(
                strsplit(dfs$ID_source, "_"), function(xx) xx[1], character(1))
            dft$prefix <- dft$ID_target
        } else if (target_platform %in% c("EPICv2", "MSA") &&
                   source_platform %in% c("EPIC", "HM450", "HM27")) {
            dfs$prefix <- dfs$ID_source
            dft$prefix <- vapply(
                strsplit(dft$ID_target, "_"), function(xx) xx[1], character(1))
        } else {
            dfs$prefix <- dfs$ID_source
            dft$prefix <- dft$ID_target
        }
        mapping <- dplyr::full_join(dfs, dft, by="prefix")
    }
    
    if (target_uniq) {
        m <- dplyr::distinct(mapping, .data[["ID_target"]], .keep_all = TRUE)
        mapping <- rbind(
            m[!is.na(m$ID_target),],
            mapping[is.na(mapping$ID_target),])
    }

    if (!include_new) { mapping <- mapping[!is.na(mapping$ID_source),] }
    if (!include_old) { mapping <- mapping[!is.na(mapping$ID_target),] }
    if (return_mapping) {
        mapping
    } else {
        stats::setNames(mapping$ID_target, mapping$ID_source)
    }
}

#' Impute of missing data of specific platform
#' 
#' @param betas named vector of beta values
#' @param platform platform
#' @param celltype celltype/tissue context of imputation, if not given, will
#' use nearest neighbor to determine.
#' @return imputed data, vector or matrix
imputeBetas <- function(betas, platform=NULL, celltype=NULL) {

    df <- sesameDataGet(sprintf("%s.imputationDefault", platform))
    d2q <- match(names(betas), df$Probe_ID)
    celltype <- names(which.max(vapply(df$data, function(x) cor(
        betas, x$median[d2q], use="na.or.complete"), numeric(1))))
    if (is.null(celltype)) {
        celltype <- "Blood"
    }
    idx <- is.na(betas)
    betas[idx] <- df$data[[celltype]]$median[d2q][idx]
    betas
}
