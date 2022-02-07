
#' Extract the probe type field from probe ID
#' This only works with the new probe ID system.
#' See https://github.com/zhou-lab/InfiniumAnnotation for illustration
#'
#' @param Probe_ID Probe ID
#' @return a vector of '1' and '2' suggesting Infinium-I and Infinium-II
#' @examples
#' probeID_designType("cg36609548_TC21")
#' @export
probeID_designType <- function(Probe_ID) {
    stopifnot(all(grepl('_', Probe_ID))) # make sure it's the new ID system
    vapply(Probe_ID, function(x) substr(
        strsplit(x,'_')[[1]][2],3,3), character(1))
}

#' Whether the probe ID is the uniq probe ID like in
#' the mouse array, e.g., cg36609548
#' 
#' @param Probe_ID Probe ID
#' @return a logical(1), whether the probe ID is based on the new ID system
isUniqProbeID <- function(Probe_ID) {
    all(grepl('_',Probe_ID))
}

#' Extract the first design category
#'
#' @param design_str Design string in e.g., the mouse array
#' @return a character vector for the design category
#' @import stringr
extractDesign <- function(design_str) {
    vapply(
        stringr::str_split(design_str, ','),
        function(x) stringr::str_split(x[[1]],';')[[1]][1], character(1))
}

#' Convert beta-value to M-value
#'
#' Logit transform a beta value vector to M-value vector.
#'
#' Convert beta-value to M-value (aka logit transform)
#' @param b vector of beta values
#' @return a vector of M values
#' @examples
#' BetaValueToMValue(c(0.1, 0.5, 0.9))
#' @export
BetaValueToMValue <- function(b) {
    log2(b/(1-b))
}

#' Convert M-value to beta-value
#'
#' Convert M-value to beta-value (aka inverse logit transform)
#' 
#' @param m a vector of M values
#' @return a vector of beta values
#' @examples
#' MValueToBetaValue(c(-3, 0, 3))
#' @export
MValueToBetaValue <- function(m) {
    2^m/(1+2^m)
}

pkgTest <- function(x) {
    if (!require(x, character.only = TRUE)) {
        stop("Optional package ", x, " not found.
Please install before continue.")
    }
}

#' Get most variable probes
#'
#' @param betas beta value matrix (row: cpg; column: sample)
#' @param n number of most variable probes
#' @return beta value matrix for the most variable probes
#' @examples
#' ## get most variable autosome probes
#' betas <- sesameDataGet('HM450.10.TCGA.PAAD.normal')
#' betas.most.variable <- bSubMostVariable(
#'     betas[names(sesameData_getAutosomeProbes('HM450')),],2000)
#' @export
bSubMostVariable <- function(betas, n=2000) {
    std <- apply(betas, 1, sd, na.rm=TRUE)
    betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
}

#' subset beta value matrix by probes
#' 
#' @param betas beta value matrix
#' @param probes probe set
#' @return subsetted beta value matrix
#' @examples
#' probes <- names(sesameData_getAutosomeProbes('HM450'))
#' betas <- sesameDataGet('HM450.1.TCGA.PAAD')$betas
#' betas <- bSubProbes(betas, probes)
#' @export
bSubProbes <- function(betas, probes) {
    if (is.null(dim(betas))) { # should also work for vector
        betas[intersect(names(betas), probes)]
    } else {
        betas[intersect(rownames(betas), probes),]
    }
}

#' subset beta value matrix by complete probes
#' 
#' @param betas beta value matrix
#' @return subsetted beta value matrix
#' @examples
#' betas <- sesameDataGet('HM450.1.TCGA.PAAD')$betas
#' betas <- bSubComplete(betas)
#' @export
bSubComplete <- function(betas) {
    if (is.null(dim(betas))) { # should also work for vector
        betas[!is.na(betas)]
    } else {
        betas[complete.cases(betas),]
    }
}

#' Annotate a data.frame using manifest
#'
#' @param df input data frame with Probe_ID as a column
#' @param probe_id the Probe_ID column name, default to "Probe_ID" or rownames
#' @param platform which array platform, guess from probe ID if not given
#' @param genome the genome build, use default if not given
#' @return a new data.frame with manifest attached
#' @examples
#' df <- data.frame(Probe_ID = c("cg00101675_BC21", "cg00116289_BC21"))
#' attachManifest(df)
#' @export
attachManifest <- function(
    df, probe_id="Probe_ID", platform=NULL, genome=NULL) {
    stopifnot(is(df, "data.frame"))
    stopifnot(probe_id %in% colnames(df))

    if (is.null(platform)) {
        platform <- inferPlatformFromProbeIDs(df[[probe_id]]) }

    genome <- sesameData_check_genome(genome, platform)

    mft <- sesameDataGet(sprintf("%s.%s.manifest", platform, genome))
    cbind(df, as.data.frame(mft)[match(df[[probe_id]], names(mft)),])
}

