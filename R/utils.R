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
    mft <- sesameDataGet(sprintf('%s.hg38.manifest', platform))
    names(mft)[as.character(seqnames(mft)) %in% chrms]
}

#' Get autosome probes
#'
#' @param platform 'EPIC', 'HM450' etc.
#' @return a vector of autosome probes
#' @examples
#' auto.probes <- getAutosomeProbes('EPIC')
#' @export
getAutosomeProbes <- function(
    platform=c('EPIC','HM450','Mouse'),
    refversion=c('hg19','hg38','mm10')) {
    
    mft <- sesameDataGet(sprintf('%s.hg38.manifest', platform))
    names(mft)[!(as.character(seqnames(mft)) %in% c('chrX', 'chrY'))]
}

#' Get most variable probes
#'
#' @param betas beta value matrix (row: cpg; column: sample)
#' @param n number of most variable probes
#' @return beta value matrix for the most variable probes
#' @examples
#' ## get most variable autosome probes
#' betas <- sesameDataGet('HM450.1.TCGA.PAAD')$betas
#' betas.most.variable <- getMostVariableProbes(
#'     betas[getAutosomeProbes('HM450'),],2000)
#' @export
getMostVariableProbes <- function(betas, n=2000) {
    sd <- apply(betas, 1, sd, na.rm=T)
    betas[names(sort(sd, decreasing=T)[1:n]),]
}

#' subset beta value matrix by probes
#' 
#' @param betas beta value matrix
#' @param probes probe set
#' @examples
#' probes <- getAutosomeProbes('HM450')
#' betas <- sesameDataGet('HM450.1.TCGA.PAAD')$betas
#' betas <- bSubProbes(betas, probes)
#' @export
bSubProbes <- function(betas, probes) {
    betas[intersect(rownames(betas), probes),]
}

#' subset beta value matrix by complete probes
#' 
#' @param betas beta value matrix
#' @examples
#' betas <- sesameDataGet('HM450.1.TCGA.PAAD')$betas
#' betas <- bSubComplete(betas)
#' @export
bSubComplete <- function(betas) {
    betas[complete.cases(betas),]
}
