
#' Detection P-value set to all zero
#'
#' @param sset a \code{SigSet}
#' @return detection p-value set to all zero
#' @param force force rerun even if result already exists
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionZero(sset)
#' 
#' @export
detectionZero <- function(sset, force=FALSE) {

    stopifnot(is(sset, "SigSet"))
    method <- "Zero"
    if (!force && method %in% names(extra(sset)$pvals)) return(sset)
    
    if (!('pvals' %in% names(extra(sset))))
        extra(sset)[['pvals']] <- list()
    
    nms <- probeNames(sset)
    extra(sset)[['pvals']][[method]] <-
        setNames(rep(0, times = length(nms)), nms)
    
    sset
}

#' Detection P-value based on ECDF of negative control
#'
#' The function takes a \code{SigSet} as input, computes detection p-value
#' using negative control probes' empirical distribution and returns a new
#' \code{SigSet} with an updated pval slot.
#'
#' @param sset a \code{SigSet}
#' @param force force rerun even if result already exists
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPnegEcdf(sset)
#' @import methods
#' @export
detectionPnegEcdf <- function(sset, force=FALSE) {

    stopifnot(is(sset, "SigSet"))
    method <- "PnegEcdf"
    if (!force && method %in% names(extra(sset)$pvals)) return(sset)

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

    ## note: no sorting here
    if (!('pvals' %in% names(extra(sset))))
        extra(sset)[['pvals']] <- list()
    
    extra(sset)[['pvals']][[method]] <- c(pIR,pIG,pII)

    sset
}

#' Detection P-value based on ECDF of out-of-band signal
#' 
#' aka pOOBAH (p-vals by Out-Of-Band Array Hybridization)
#'
#' The function takes a \code{SigSet} as input, computes detection p-value
#' using out-of-band probes empirical distribution and returns a new
#' \code{SigSet} with an updated pval slot.
#'
#' @name detectionPoobEcdf
#' @param sset a \code{SigSet}
#' @param force force rerun even if result already exists
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPoobEcdf(sset)
#' 
#' @export
detectionPoobEcdf <- function(sset, force=FALSE) {

    stopifnot(is(sset, "SigSet"))
    method <- "pOOBAH"
    if (!force && method %in% names(extra(sset)$pvals)) return(sset)

    funcG <- ecdf(oobG(sset))
    funcR <- ecdf(oobR(sset))

    ## p-value is the minimium detection p-value of the 2 alleles
    pIR <- 1-apply(cbind(funcR(IR(sset)[,'M']), funcR(IR(sset)[,'U'])),1,max)
    pIG <- 1-apply(cbind(funcG(IG(sset)[,'M']), funcG(IG(sset)[,'U'])),1,max)
    pII <- 1-apply(cbind(funcG(II(sset)[,'M']), funcR(II(sset)[,'U'])),1,max)

    names(pIR) <- rownames(IR(sset))
    names(pIG) <- rownames(IG(sset))
    names(pII) <- rownames(II(sset))

    ## should have PoobEcdf aliased
    if (!('pvals' %in% names(extra(sset))))
        extra(sset)[['pvals']] <- list()
    extra(sset)[['pvals']][[method]] <- c(pIR,pIG,pII)
    
    sset
}

#' @rdname detectionPoobEcdf
#' @export
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- pOOBAH(sset)
pOOBAH <- detectionPoobEcdf

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
#' @param force force rerun even if result already exists
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPnegNorm(sset)
#' 
#' @export
detectionPnegNorm <- function(sset, force=FALSE) {

    stopifnot(is(sset, "SigSet"))
    method <- "PnegNorm"
    if (!force && method %in% names(extra(sset)$pvals)) return(sset)
    
    negctls <- negControls(sset)
    muG <- median(negctls$G)
    sdG <- sd(negctls$G)
    muR <- median(negctls$R)
    sdR <- sd(negctls$R)
    sset <- detectionPfixedNorm(sset, muG, sdG, muR, sdR)

    if (!('pvals' %in% names(extra(sset))))
        extra(sset)[['pvals']] <- list()
    
    extra(sset)[['pvals']][[method]] <- extra(sset)[['pvals']][['PfixNorm']]
    sset
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
#' @param force force rerun even if result already exists
#' @param sdG SD of background in Grn channel
#' @param muG mean of background in Grn channel
#' @param sdR SD of background in Red channel
#' @param muR mean of background in Red channel
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPfixedNorm(sset)
#' 
#' @export
detectionPfixedNorm <- function(
    sset, muG = 500, sdG = 200, muR = 500, sdR = 200,
    force = FALSE) {
    
    stopifnot(is(sset, "SigSet"))
    method <- "PfixedNorm"
    if (!force && method %in% names(extra(sset)$pvals)) return(sset)
    
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
    
    if (!('pvals' %in% names(extra(sset))))
        extra(sset)[['pvals']] <- list()
    
    extra(sset)[['pvals']][[method]] <- c(pIR,pIG,pII)
    sset
}

#' Detection P-value emulating Genome Studio
#'
#' The function takes a \code{SigSet} as input, computes detection p-value
#' using negative control probes parametrized in a normal distribution a la
#' Genome Studio and returns a new \code{SigSet} with an updated pval slot.
#' 
#' P-value is calculated using negative control probes as the estimate of
#' background where Grn channel and Red channel are merged. But when
#' estimating p-value the Red and Grn are summed (non-ideal).
#'
#' @param sset a \code{SigSet}
#' @param force force rerun even if result already exists
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPnegNormGS(sset)
#' 
#' @export
detectionPnegNormGS <- function(sset, force=FALSE) {
    
    stopifnot(is(sset, "SigSet"))
    method <- "PnegNormGS"
    if (!force && method %in% names(extra(sset)$pvals)) return(sset)
    
    negctls <- negControls(sset)
    BGsd <- sd(c(negctls$G, negctls$R))
    BGmu <- mean(c(negctls$G, negctls$R))
    pIR <- 1 - pnorm(rowSums(IR(sset)), mean=BGmu, sd=BGsd)
    pIG <- 1 - pnorm(rowSums(IG(sset)), mean=BGmu, sd=BGsd)
    pII <- 1 - pnorm(rowSums(II(sset)), mean=BGmu, sd=BGsd)
    names(pIR) <- rownames(IR(sset))
    names(pIG) <- rownames(IG(sset))
    names(pII) <- rownames(II(sset))
    
    if (!('pvals' %in% names(extra(sset))))
        extra(sset)[['pvals']] <- list()
    
    extra(sset)[['pvals']][[method]] <- c(pIR,pIG,pII)
    sset
    
}

#' Detection P-value based on normal fitting the negative controls,
#' channels are first summed
#'
#' The function takes a \code{SigSet} as input, computes detection p-value
#' using negative control probes parametrized in a normal distribution with
#' the two channels summed first and returns a new \code{SigSet} with an
#' updated pval slot. The SD is summed to emulate the SD of the summed
#' signal (not the most accurate treatment).
#'
#' @param sset a \code{SigSet}
#' @param force force rerun even if result already exists
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPnegNormTotal(sset)
#' 
#' @export
detectionPnegNormTotal <- function(sset, force=FALSE) {
    
    ## sort of how minfi does it (sd instead of MAD)
    stopifnot(is(sset, "SigSet"))
    method <- "PnegNormTotal"
    if (!force && method %in% names(extra(sset)$pvals)) return(sset)
    
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
    
    if (!('pvals' %in% names(extra(sset))))
        extra(sset)[['pvals']] <- list()
    
    extra(sset)[['pvals']][[method]] <- c(pIR,pIG,pII)
    sset
}

