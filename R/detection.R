
#' Detection P-value based on ECDF of negative control
#'
#' The function takes a \code{SigSet} as input, computes detection p-value
#' using negative control probes' empirical distribution and returns a new
#' \code{SigSet} with an updated mask slot.
#'
#' @param sset a \code{SigSet}
#' @param pval.threshold minimum p-value to mask
#' @param return.pval whether to return p-values, instead of a
#' masked \code{SigSet}
#' @return a \code{SigSet}, or a p-value vector if return.pval is TRUE
#' @examples
#' sset <- sesameDataGet("HM450.1.TCGA.PAAD")$sset
#' sum(mask(sset))
#' sset_with_mask <- detectionPnegEcdf(sset)
#' sum(mask(sset_with_mask))
#' @import methods
#' @export
detectionPnegEcdf <- function(sset, return.pval = FALSE, pval.threshold=0.05) {

    stopifnot(is(sset, "SigSet"))

    negctls <- negControls(sset)
    funcG <- ecdf(negctls$G)
    funcR <- ecdf(negctls$R)

    ## p-value is the minimium detection p-value of the 2 alleles
    pIR <- 1-apply(cbind(funcR(IR(sset)[,'M']), funcR(IR(sset)[,'U'])),1,max)
    pIG <- 1-apply(cbind(funcG(IG(sset)[,'M']), funcG(IG(sset)[,'U'])),1,max)
    pII <- 1-apply(cbind(funcG(II(sset)[,'M']), funcR(II(sset)[,'U'])),1,max)

    names(pIR) <- rownames(IR(sset))
    names(pIG) <- rownames(IG(sset))
    names(pII) <- rownames(II(sset))

    pvals <- c(pIR,pIG,pII)
    if (return.pval) {
        return(pvals[order(names(pvals))])
    }

    addMask(sset, pvals > pval.threshold)
}

#' Detection P-value based on ECDF of out-of-band signal
#' 
#' aka pOOBAH (p-vals by Out-Of-Band Array Hybridization)
#'
#' The function takes a \code{SigSet} as input, computes detection p-value
#' using out-of-band probes empirical distribution and returns a new
#' \code{SigSet} with an updated mask slot.
#'
#' @name detectionPoobEcdf
#' @param sset a \code{SigSet}
#' @param pval.threshold minimum p-value to mask
#' @param return.pval whether to return p-values, instead of a
#' masked \code{SigSet}
#' @return a \code{SigSet}, or a p-value vector if return.pval is TRUE
#' @examples
#' sset <- sesameDataGet("HM450.1.TCGA.PAAD")$sset
#' sum(mask(sset))
#' sset_with_mask <- detectionPoobEcdf(sset)
#' sum(mask(sset_with_mask))
#' @export
detectionPoobEcdf <- function(sset, return.pval = FALSE, pval.threshold=0.05) {

    stopifnot(is(sset, "SigSet"))

    funcG <- ecdf(oobG(sset))
    funcR <- ecdf(oobR(sset))

    ## p-value is the minimium detection p-value of the 2 alleles
    pIR <- 1-apply(cbind(funcR(IR(sset)[,'M']), funcR(IR(sset)[,'U'])),1,max)
    pIG <- 1-apply(cbind(funcG(IG(sset)[,'M']), funcG(IG(sset)[,'U'])),1,max)
    pII <- 1-apply(cbind(funcG(II(sset)[,'M']), funcR(II(sset)[,'U'])),1,max)

    names(pIR) <- rownames(IR(sset))
    names(pIG) <- rownames(IG(sset))
    names(pII) <- rownames(II(sset))

    pvals <- c(pIR,pIG,pII)
    if (return.pval) {
        return(pvals[order(names(pvals))])
    }

    addMask(sset, pvals > pval.threshold)
}

#' Detection P-value based on ECDF of out-of-band signal
#' 
#' aka pOOBAH2 (p-vals by Out-Of-Band Array Hybridization)
#'
#' The function takes a \code{SigSet} as input, computes detection p-value
#' using out-of-band probes empirical distribution and returns a new
#' \code{SigSet} with an updated mask slot.
#'
#' The difference between this function and the original pOOBAH
#' is that pOOBAH2 is based on background-subtracted and dyebias
#' corrected signal and do not distinguish the color channel difference.
#'
#' @name detectionPoobEcdf2
#' @param sset a \code{SigSet}
#' @param pval.threshold minimum p-value to mask
#' @param return.pval whether to return p-values, instead of a
#' masked \code{SigSet}
#' @return a \code{SigSet}, or a p-value vector if return.pval is TRUE
#' @examples
#' sset <- sesameDataGet("HM450.1.TCGA.PAAD")$sset
#' sum(mask(sset))
#' sset_with_mask <- detectionPoobEcdf(sset)
#' sum(mask(sset_with_mask))
#' @export
detectionPoobEcdf2 <- function(sset, return.pval = FALSE, pval.threshold=0.05){

    stopifnot(is(sset, "SigSet"))

    func <- ecdf(c(oobG(sset), oobR(sset)))

    ## p-value is the minimium detection p-value of the 2 alleles
    pIR <- 1-apply(cbind(func(IR(sset)[,'M']), func(IR(sset)[,'U'])),1,max)
    pIG <- 1-apply(cbind(func(IG(sset)[,'M']), func(IG(sset)[,'U'])),1,max)
    pII <- 1-apply(cbind(func(II(sset)[,'M']), func(II(sset)[,'U'])),1,max)

    names(pIR) <- rownames(IR(sset))
    names(pIG) <- rownames(IG(sset))
    names(pII) <- rownames(II(sset))

    pvals <- c(pIR,pIG,pII)
    if (return.pval) {
        return(pvals[order(names(pvals))])
    }

    addMask(sset, pvals > pval.threshold)
}

#' @rdname detectionPoobEcdf
#' @export
#' @examples
#' sset <- sesameDataGet("HM450.1.TCGA.PAAD")$sset
#' sum(mask(sset))
#' sset_with_mask <- pOOBAH(sset)
#' sum(mask(sset_with_mask))
pOOBAH <- detectionPoobEcdf

#' @rdname detectionPoobEcdf2
#' @export
#' @examples
#' sset <- sesameDataGet("HM450.1.TCGA.PAAD")$sset
#' sum(mask(sset))
#' sset_with_mask <- pOOBAH2(sset)
#' sum(mask(sset_with_mask))
pOOBAH2 <- detectionPoobEcdf2

#################################################
## Detection P-value based on normal distribution
#################################################

#' Detection P-value based on normal fitting the negative controls
#'
#' The function takes a \code{SigSet} as input, computes detection p-value
#' using negative control probes parametrized in a normal distribution and
#' returns a new \code{SigSet} with an updated pval slot.
#' 
#' Background of Grn and Red are estimated separately from negative control
#' probes-parameterized normal distribution. p-value is taken from the
#' minimum of the p-value of the two alleles (color depends on probe design).
#'
#' @param sset a \code{SigSet}
#' @param pval.threshold minimum p-value to mask
#' @param return.pval whether to return p-values, instead of a
#' masked \code{SigSet}
#' @return a \code{SigSet}, or a p-value vector if return.pval is TRUE
#' @examples
#' sset <- sesameDataGet("HM450.1.TCGA.PAAD")$sset
#' sum(mask(sset))
#' sset_with_mask <- detectionPnegNorm(sset)
#' sum(mask(sset_with_mask))
#' @export
detectionPnegNorm <- function(sset,
    pval.threshold=0.05, return.pval = FALSE) {

    stopifnot(is(sset, "SigSet"))
    
    negctls <- negControls(sset)
    muG <- median(negctls$G)
    sdG <- sd(negctls$G)
    muR <- median(negctls$R)
    sdR <- sd(negctls$R)
    detectionPfixedNorm(sset, muG, sdG, muR, sdR,
        return.pval = return.pval, pval.threshold = 0.05)
}

#' Detection P-value based on normal fitting with gived parameters
#'
#' The function takes a \code{SigSet} as input, computes detection p-value
#' using negative control probes parametrized in a normal distribution and
#' returns a new \code{SigSet} with an updated pval slot.
#' 
#' Background of Grn and Red are estimated separately from a fixed normal 
#' distribution. p-value is taken from the minimum of the p-value of the 
#' two alleles (color depends on probe design).
#'
#' @param sset a \code{SigSet}
#' @param sdG SD of background in Grn channel
#' @param muG mean of background in Grn channel
#' @param sdR SD of background in Red channel
#' @param muR mean of background in Red channel
#' @param pval.threshold minimum p-value to mask
#' @param return.pval whether to return p-values, instead of a
#' masked \code{SigSet}
#' @return detection p-value
#' @examples
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' sum(mask(sset))
#' sset_with_mask <- detectionPfixedNorm(sset)
#' sum(mask(sset_with_mask))
#' 
#' @export
detectionPfixedNorm <- function(
    sset, muG = 500, sdG = 200, muR = 500, sdR = 200,
    pval.threshold=0.05, return.pval = FALSE) {
    
    stopifnot(is(sset, "SigSet"))
    
    pIR <- 1 - pnorm(
        pmax(IR(sset)[,'M'], IR(sset)[,'U']), mean=muR, sd=sdR)
    pIG <- 1 - pnorm(
        pmax(IG(sset)[,'M'], IG(sset)[,'U']), mean=muG, sd=sdG)
    pII <- pmin(
        1 - pnorm(II(sset)[,'M'], mean=muG, sd=sdG),
        1 - pnorm(II(sset)[,'U'], mean=muR, sd=sdR))
    
    names(pIR) <- rownames(IR(sset))
    names(pIG) <- rownames(IG(sset))
    names(pII) <- rownames(II(sset))
    
    pvals <- c(pIR,pIG,pII)
    if (return.pval) {
        return(pvals[order(names(pvals))])
    }

    addMask(sset, pvals > pval.threshold)
}

#' Detection P-value emulating Genome Studio
#'
#' The function takes a \code{SigSet} as input, computes detection p-value
#' using negative control probes parametrized in a normal distribution a la
#' Genome Studio and returns a new \code{SigSet} with an updated mask slot.
#' 
#' P-value is calculated using negative control probes as the estimate of
#' background where Grn channel and Red channel are merged. But when
#' estimating p-value the Red and Grn are summed (non-ideal).
#'
#' @param sset a \code{SigSet}
#' @param pval.threshold minimum p-value to mask
#' @param return.pval whether to return p-values, instead of a
#' masked \code{SigSet}
#' @return detection p-value
#' @examples
#' sset <- sesameDataGet("HM450.1.TCGA.PAAD")$sset
#' sum(mask(sset))
#' sset_with_mask <- detectionPnegNormGS(sset)
#' sum(mask(sset_with_mask))
#' @export
detectionPnegNormGS <- function(sset,
    pval.threshold=0.05, return.pval = FALSE) {
    
    stopifnot(is(sset, "SigSet"))
    
    negctls <- negControls(sset)
    BGsd <- sd(c(negctls$G, negctls$R))
    BGmu <- mean(c(negctls$G, negctls$R))
    pIR <- 1 - pnorm(rowSums(IR(sset)), mean=BGmu, sd=BGsd)
    pIG <- 1 - pnorm(rowSums(IG(sset)), mean=BGmu, sd=BGsd)
    pII <- 1 - pnorm(rowSums(II(sset)), mean=BGmu, sd=BGsd)
    names(pIR) <- rownames(IR(sset))
    names(pIG) <- rownames(IG(sset))
    names(pII) <- rownames(II(sset))
    
    pvals <- c(pIR,pIG,pII)
    if (return.pval) {
        return(pvals[order(names(pvals))])
    }

    addMask(sset, pvals > pval.threshold)
}

#' Detection P-value based on normal fitting the negative controls,
#' channels are first summed
#'
#' The function takes a \code{SigSet} as input, computes detection p-value
#' using negative control probes parametrized in a normal distribution with
#' the two channels summed first and returns a new \code{SigSet} with an
#' updated mask slot. The SD is summed to emulate the SD of the summed
#' signal (not the most accurate treatment).
#'
#' @param sset a \code{SigSet}
#' @param pval.threshold minimum p-value to mask
#' @param return.pval whether to return p-values, instead of a
#' masked \code{SigSet}
#' @return detection p-value
#' @examples
#' sset <- sesameDataGet("HM450.1.TCGA.PAAD")$sset
#' sum(mask(sset))
#' sset_with_mask <- detectionPnegNormTotal(sset)
#' sum(mask(sset_with_mask))
#' @export
detectionPnegNormTotal <- function(sset,
    pval.threshold=0.05,  return.pval = FALSE) {
    
    ## sort of how minfi does it (sd instead of MAD)
    stopifnot(is(sset, "SigSet"))
    
    negctls <- negControls(sset)
    sdG <- sd(negctls$G)
    muG <- median(negctls$G)
    sdR <- sd(negctls$R)
    muR <- median(negctls$R)
    pIR <- 1 - pnorm(rowSums(IR(sset)), mean=2*muR, sd=2*sdR)
    pIG <- 1 - pnorm(rowSums(IG(sset)), mean=2*muG, sd=2*sdG)
    pII <- 1 - pnorm(rowSums(II(sset)), mean=muR+muG, sd=sdR+sdG)
    names(pIR) <- rownames(IR(sset))
    names(pIG) <- rownames(IG(sset))
    names(pII) <- rownames(II(sset))
    
    pvals <- c(pIR,pIG,pII)
    if (return.pval) {
        return(pvals[order(names(pvals))])
    }

    addMask(sset, pvals > pval.threshold)
}

