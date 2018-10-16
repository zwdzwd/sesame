#' Detection P-value based on ECDF of negative control
#'
#' The function takes a \code{SigSet} as input, computes detection p-value
#' using negative control probes' empirical distribution and returns a new
#' \code{SigSet} with an updated pval slot.
#'
#' @param sset a \code{SigSet}
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPnegEcdf(sset)
#' 
#' @export
detectionPnegEcdf <- function(sset) {

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

    pval <- c(pIR,pIG,pII)
    pval(sset) <- pval[order(names(pval))]
    sset
    
}

#' Detection P-value based on normal fitting the negative controls
#'
#' The function takes a \code{SigSet} as input, computes detection p-value
#' using negative control probes parametrized in a normal distribution and
#' returns a new \code{SigSet} with an updated pval slot.
#'
#' @param sset a \code{SigSet}
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPnegNorm(sset)
#' 
#' @export
detectionPnegNorm <- function(sset) {

    stopifnot(is(sset, "SigSet"))
    negctls <- negControls(sset)
    sdG <- sd(negctls$G)
    muG <- median(negctls$G)
    sdR <- sd(negctls$R)
    muR <- median(negctls$R)
    pIR <- 1 - pnorm(pmax(IR(sset)[,'M'], IR(sset)[,'U']), mean=muR, sd=sdR)
    pIG <- 1 - pnorm(pmax(IG(sset)[,'M'], IG(sset)[,'U']), mean=muG, sd=sdG)
    pII <- pmin(
        1 - pnorm(II(sset)[,'M'], mean=muG, sd=sdG),
        1 - pnorm(II(sset)[,'U'], mean=muR, sd=sdR))

    names(pIR) <- rownames(IR(sset))
    names(pIG) <- rownames(IG(sset))
    names(pII) <- rownames(II(sset))

    pval <- c(pIR,pIG,pII)
    pval(sset) <- pval[order(names(pval))]
    sset
    
}

#' Detection P-value emulating Genome Studio
#'
#' The function takes a \code{SigSet} as input, computes detection p-value
#' using negative control probes parametrized in a normal distribution a la
#' Genome Studio and returns a new \code{SigSet} with an updated pval slot.
#'
#' @param sset a \code{SigSet}
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPnegNormGS(sset)
#' 
#' @export
detectionPnegNormGS <- function(sset) {

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
    pval <- c(pIR,pIG,pII)
    pval(sset) <- pval[order(names(pval))]
    sset
    
}

#' Detection P-value based on normal fitting the negative controls,
#' channels are first summed
#'
#' The function takes a \code{SigSet} as input, computes detection p-value
#' using negative control probes parametrized in a normal distribution with
#' the two channels summed first and returns a new \code{SigSet} with an
#' updated pval slot.
#'
#' @param sset a \code{SigSet}
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPnegNormTotal(sset)
#' 
#' @export
detectionPnegNormTotal <- function(sset) {

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

    pval <- c(pIR,pIG,pII)
    pval(sset) <- pval[order(names(pval))]
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
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPoobEcdf(sset)
#' 
#' @export
detectionPoobEcdf <- function(sset) {

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

    pval <- c(pIR,pIG,pII)
    pval(sset) <- pval[order(names(pval))]
    sset
}

#' @rdname detectionPoobEcdf
pOOBAH <- detectionPoobEcdf
