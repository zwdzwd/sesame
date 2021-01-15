
#' Extract the probe type field from probe ID
#' This only works with the new probe ID system.
#' See https://github.com/zhou-lab/InfiniumAnnotation for illustration
#'
#' @param Probe_ID Probe ID
#' @return a vector of '1' and '2' suggesting Infinium-I and Infinium-II
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

## print message
smessage <- function(...) {
    cat('[', as.character(Sys.time()),'] ', ..., '\n', sep='')
}

pkgTest <- function(x) {
    if (!require(x, character.only = TRUE)) {
        stop("Optional package ", x, " not found.
Please install before continue.")
    }
}

#' Get probes by genomic region
#'
#' The function takes a genomic coordinate and output the a vector of probes
#' on the specified platform that falls in the given genomic region.
#'
#' @param chrm chromosome
#' @param beg begin, 1 if omitted
#' @param end end, chromosome end if omitted
#' @param platform EPIC or HM450
#' @param refversion hg38 or hg19
#' @return probes that fall into the given region
#' @importMethodsFrom IRanges subsetByOverlaps
#' @examples
#' getProbesByRegion('chr5', 135413937, 135419936,
#'     refversion = 'hg19', platform = 'HM450')
#' @export
getProbesByRegion <- function(
    chrm, beg = 1, end = -1,
    platform = c('EPIC','HM450'),
    refversion = c('hg38','hg19')) {

    platform <- match.arg(platform)
    refversion <- match.arg(refversion)
    
    if (end < 0) {
        end <- sesameDataGet(paste0(
            'genomeInfo.', refversion))$seqInfo[chrm]@seqlengths
    }

    probes <- sesameDataGet(paste0(
        platform, '.probeInfo'))[[paste0('mapped.probes.', refversion)]]
    
    if (!(chrm %in% GenomicRanges::seqinfo(probes)@seqnames)) {
        stop('No probes found in this reference');
    }
    message(sprintf('Extracting probes from %s:%d-%d.\n', chrm, beg, end))
    target.region <- GenomicRanges::GRanges(chrm, IRanges::IRanges(beg, end))
    subsetByOverlaps(probes, target.region)
}

#' Get Probes by Chromosome
#'
#' @param chrms chromosomes to subset
#' @param platform EPIC, HM450, Mouse
#' @param refversion hg19, hg38, mm10
#' @return a vector of probes on the selected chromosomes
#' @examples
#' sex.probes <- getProbesByChromosome(c('chrX','chrY'))
#' @export
getProbesByChromosome <- function(
    chrms, platform = c('EPIC','HM450'),
    refversion=c('hg19','hg38')) {

    platform <- match.arg(platform)
    refversion <- match.arg(refversion)
    
    mft <- sesameDataGet(sprintf('%s.%s.manifest', platform, refversion))
    names(mft)[as.character(GenomicRanges::seqnames(mft)) %in% chrms]
}

#' Get autosome probes
#'
#' @param platform 'EPIC', 'HM450' etc.
#' @param refversion hg19, hg38, mm10
#' @return a vector of autosome probes
#' @examples
#' auto.probes <- getAutosomeProbes('EPIC')
#' @export
getAutosomeProbes <- function(
    platform=c('EPIC','HM450','Mouse'),
    refversion=c('hg19','hg38','mm10')) {

    platform <- match.arg(platform)
    refversion <- match.arg(refversion)
    
    mft <- sesameDataGet(sprintf(
        '%s.%s.manifest', platform, refversion))
    names(mft)[!(as.character(GenomicRanges::seqnames(mft)) %in% c('chrX', 'chrY'))]
}

#' Get most variable probes
#'
#' @param betas beta value matrix (row: cpg; column: sample)
#' @param n number of most variable probes
#' @return beta value matrix for the most variable probes
#' @examples
#' ## get most variable autosome probes
#' betas <- sesameDataGet('HM450.10.TCGA.PAAD.normal')
#' betas.most.variable <- getMostVariableProbes(
#'     betas[getAutosomeProbes('HM450'),],2000)
#' @export
getMostVariableProbes <- function(betas, n=2000) {
    std <- apply(betas, 1, sd, na.rm=TRUE)
    betas[names(sort(std, decreasing=TRUE)[seq_len(n)]),]
}

#' subset beta value matrix by probes
#' 
#' @param betas beta value matrix
#' @param probes probe set
#' @return subsetted beta value matrix
#' @examples
#' probes <- getAutosomeProbes('HM450')
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

extraHas <- function(sset, k) {
    k %in% names(extra(sset))
}

extraGet <- function(sset, k) {
    extra(sset)[[k]]
}

extraSet <- function(sset, k, v) {
    extra(sset)[[k]] <- v
    sset
}

