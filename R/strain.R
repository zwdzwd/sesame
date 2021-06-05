
#' Infer strain information for mouse array
#'
#' @param vafs Variant allele frequency vector
#' @param strain_snp_table if not given download the default from sesameData
#' @return a list of best guess, p-value of the best guess
#' and the probabilities of all strains
#' @examples
#' sesameDataCache("MM285") # if not done yet
#' sdf <- sesameDataGet('MM285.1.SigDF')
#' vafs <- betaToAF(getBetas(dyeBiasNL(noob(sdf))))
#' inferStrain(vafs)
#' @import tibble
#' @export
inferStrain <- function(vafs, strain_snp_table = NULL) {

    if (is.null(strain_snp_table)) {
        strain_snp_table <- sesameDataGet('MM285.strain.snp.table')
    }
    vafs[is.na(vafs)] <- 0.5
    probes <- intersect(names(vafs), rownames(strain_snp_table))
    bb <- vapply(probes, function(p) {
        dnorm(vafs[p], mean=strain_snp_table[p,], sd=0.8)
    }, numeric(ncol(strain_snp_table)))
    bbloglik <- apply(bb,1,function(x) sum(log(x),na.rm=TRUE))
    probs <- exp(bbloglik - max(bbloglik))
    best.index <- which.max(probs)

    list(
        best = names(best.index),
        pval = (sum(probs) - probs[best.index]) / sum(probs),
        probs = probs/sum(probs))
}

#' convert betas to variant allele frequency
#'
#' @param betas beta value
#' @return SNP variant allele frequency
#' @examples
#' sesameDataCache("MM285") # if not done yet
#' sdf <- sesameDataGet('MM285.1.SigDF')
#' vafs <- betaToAF(getBetas(dyeBiasNL(noob(sdf))))
#' @export
betaToAF <- function(betas) {

    mft <- sesameDataGet('MM285.mm10.manifest')
    mft <- mft[grep('^rs', names(mft))]
    design <- GenomicRanges::mcols(mft)[['design']]
    toFlip <- !setNames(as.logical(substr(
        design, nchar(design), nchar(design))), names(mft))

    vaf <- betas[grep('rs', names(betas))]
    vaf[toFlip[names(vaf)]] <- 1-vaf[toFlip[names(vaf)]]
    vaf
}
