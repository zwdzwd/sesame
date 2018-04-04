#' Detection P-value based on ECDF of negative control
#'
#' Returns a new \code{SigSet} with updated pval slot
#'
#' @param sset a \code{SigSet}
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPnegEcdf(sset)
#' 
#' @export
detectionPnegEcdf <- function(sset) {
    
    negctls <- negControls(sset)
    funcG <- ecdf(negctls$G)
    funcR <- ecdf(negctls$R)

    ## p-value is the minimium detection p-value of the 2 alleles
    pIR <- 1-apply(cbind(funcR(sset@IR[,'M']), funcR(sset@IR[,'U'])),1,max)
    pIG <- 1-apply(cbind(funcG(sset@IG[,'M']), funcG(sset@IG[,'U'])),1,max)
    pII <- 1-apply(cbind(funcG(sset@II[,'M']), funcR(sset@II[,'U'])),1,max)

    names(pIR) <- rownames(sset@IR)
    names(pIG) <- rownames(sset@IG)
    names(pII) <- rownames(sset@II)

    pval <- c(pIR,pIG,pII)
    sset@pval <- pval[order(names(pval))]
    sset
    
}

#' Detection P-value based on normal fitting the negative controls
#'
#' Returns a new \code{SigSet} with updated pval slot
#'
#' @param sset a \code{SigSet}
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPnegNorm(sset)
#' 
#' @export
detectionPnegNorm <- function(sset) {
    
    negctls <- negControls(sset)
    sdG <- sd(negctls$G)
    muG <- median(negctls$G)
    sdR <- sd(negctls$R)
    muR <- median(negctls$R)
    pIR <- 1 - pnorm(pmax(sset@IR[,'M'], sset@IR[,'U']), mean=muR, sd=sdR)
    pIG <- 1 - pnorm(pmax(sset@IG[,'M'], sset@IG[,'U']), mean=muG, sd=sdG)
    pII <- pmin(
        1 - pnorm(sset@II[,'M'], mean=muG, sd=sdG),
        1 - pnorm(sset@II[,'U'], mean=muR, sd=sdR))

    names(pIR) <- rownames(sset@IR)
    names(pIG) <- rownames(sset@IG)
    names(pII) <- rownames(sset@II)

    pval <- c(pIR,pIG,pII)
    sset@pval <- pval[order(names(pval))]
    sset
    
}

#' Detection P-value emulating Genome Studio
#'
#' Returns a new \code{SigSet} with updated pval slot
#'
#' @param sset a \code{SigSet}
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPnegNormGS(sset)
#' 
#' @export
detectionPnegNormGS <- function(sset) {
    
    negctls <- negControls(sset)
    BGsd <- sd(c(negctls$G, negctls$R))
    BGmu <- mean(c(negctls$G, negctls$R))
    pIR <- 1 - pnorm(rowSums(sset@IR), mean=BGmu, sd=BGsd)
    pIG <- 1 - pnorm(rowSums(sset@IG), mean=BGmu, sd=BGsd)
    pII <- 1 - pnorm(rowSums(sset@II), mean=BGmu, sd=BGsd)
    names(pIR) <- rownames(sset@IR)
    names(pIG) <- rownames(sset@IG)
    names(pII) <- rownames(sset@II)
    pval <- c(pIR,pIG,pII)
    sset@pval <- pval[order(names(pval))]
    sset
    
}

#' Detection P-value based on normal fitting the negative controls,
#' channels are first summed
#'
#' Returns a new \code{SigSet} with updated pval slot
#'
#' @param sset a \code{SigSet}
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPnegNormTotal(sset)
#' 
#' @export
detectionPnegNormTotal <- function(sset) {
    
    ## how minfi does it
    negctls <- negControls(sset)
    sdG <- sd(negctls$G)
    muG <- median(negctls$G)
    sdR <- sd(negctls$R)
    muR <- median(negctls$R)
    pIR <- 1 - pnorm(rowSums(sset@IR), mean=2*muR, sd=2*sdR)
    pIG <- 1 - pnorm(rowSums(sset@IG), mean=2*muG, sd=2*sdG)
    pII <- 1 - pnorm(rowSums(sset@II), mean=muR+muG, sd=sdR+sdG)
    names(pIR) <- rownames(sset@IR)
    names(pIG) <- rownames(sset@IG)
    names(pII) <- rownames(sset@II)

    pval <- c(pIR,pIG,pII)
    sset@pval <- pval[order(names(pval))]
    sset
}

#' Detection P-value based on ECDF of out-of-band signal
#'
#' Returns a new \code{SigSet} with updated pval slot
#'
#' @param sset a \code{SigSet}
#' @return detection p-value
#' @examples
#' sset <- makeExampleSeSAMeDataSet()
#' sset <- detectionPoobEcdf(sset)
#' 
#' @export
detectionPoobEcdf <- function(sset) {
    
    funcG <- ecdf(sset@oobG)
    funcR <- ecdf(sset@oobR)

    ## p-value is the minimium detection p-value of the 2 alleles
    pIR <- 1-apply(cbind(funcR(sset@IR[,'M']), funcR(sset@IR[,'U'])),1,max)
    pIG <- 1-apply(cbind(funcG(sset@IG[,'M']), funcG(sset@IG[,'U'])),1,max)
    pII <- 1-apply(cbind(funcG(sset@II[,'M']), funcR(sset@II[,'U'])),1,max)

    names(pIR) <- rownames(sset@IR)
    names(pIG) <- rownames(sset@IG)
    names(pII) <- rownames(sset@II)

    pval <- c(pIR,pIG,pII)
    sset@pval <- pval[order(names(pval))]
    sset
}

