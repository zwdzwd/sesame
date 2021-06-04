## get negative control probes
negControls <- function(sdf) {
    stopifnot(is(sdf, "SigDF"))
    negctls <- controls(sdf)[grep(
        'negative', tolower(rownames(controls(sdf)))),]
    negctls <- subset(negctls, col!=-99)
    negctls
}

#' Detection P-value based on ECDF of negative control
#'
#' The function takes a \code{SigDF} as input, computes detection p-value
#' using negative control probes' empirical distribution and returns a new
#' \code{SigDF} with an updated mask slot.
#'
#' @param sdf a \code{SigDF}
#' @param pval.threshold minimum p-value to mask
#' @param return.pval whether to return p-values, instead of a
#' masked \code{SigDF}
#' @return a \code{SigDF}, or a p-value vector if return.pval is TRUE
#' @examples
#' sdf <- sesameDataGet("HM450.1.TCGA.PAAD")$sdf
#' sum(mask(sdf))
#' sdf_with_mask <- detectionPnegEcdf(sdf)
#' sum(mask(sdf_with_mask))
#' @import methods
#' @export
detectionPnegEcdf <- function(sdf, return.pval = FALSE, pval.threshold=0.05) {

    stopifnot(is(sdf, "SigDF"))

    negctls <- negControls(sdf)
    funcG <- ecdf(negctls$G)
    funcR <- ecdf(negctls$R)

    ## p-value is the minimium detection p-value of the 2 alleles
    pvals = setNames(pmin(
        with(sdf, 1-funcR(pmax(MR, UR, na.rm=TRUE))),
        with(sdf, 1-funcG(pmax(MG, UG, na.rm=TRUE)))), sdf$Probe_ID)
        
    if (return.pval) { return(pvals) }

    addMask(sdf, pvals > pval.threshold)
}

#' Detection P-value based on ECDF of out-of-band signal
#' 
#' aka pOOBAH (p-vals by Out-Of-Band Array Hybridization)
#'
#' The function takes a \code{SigDF} as input, computes detection p-value
#' using out-of-band probes empirical distribution and returns a new
#' \code{SigDF} with an updated mask slot.
#'
#' @name detectionPoobEcdf
#' @param sdf a \code{SigDF}
#' @param pval.threshold minimum p-value to mask
#' @param return.pval whether to return p-values, instead of a
#' masked \code{SigDF}
#' @return a \code{SigDF}, or a p-value vector if return.pval is TRUE
#' @examples
#' sdf <- sesameDataGet("HM450.1.TCGA.PAAD")$sdf
#' sum(mask(sdf))
#' sdf_with_mask <- detectionPoobEcdf(sdf)
#' sum(mask(sdf_with_mask))
#' @export
detectionPoobEcdf <- function(sdf, return.pval = FALSE, pval.threshold=0.05) {

    stopifnot(is(sdf, "SigDF"))

    funcG <- with(IR(sdf), ecdf(c(MG,UG)))
    funcR <- with(IG(sdf), ecdf(c(MR,UR)))

    ## p-value is the minimium detection p-value of the 2 alleles
    pvals = setNames(pmin(
        with(sdf, 1-funcR(pmax(MR, UR, na.rm=TRUE))),
        with(sdf, 1-funcG(pmax(MG, UG, na.rm=TRUE)))), sdf$Probe_ID)
        
    if (return.pval) { return(pvals) }

    addMask(sdf, pvals > pval.threshold)
}

#' Detection P-value based on ECDF of out-of-band signal
#' 
#' aka pOOBAH2 (p-vals by Out-Of-Band Array Hybridization)
#'
#' The function takes a \code{SigDF} as input, computes detection p-value
#' using out-of-band probes empirical distribution and returns a new
#' \code{SigDF} with an updated mask slot.
#'
#' The difference between this function and the original pOOBAH
#' is that pOOBAH2 is based on background-subtracted and dyebias
#' corrected signal and do not distinguish the color channel difference.
#'
#' @name detectionPoobEcdf2
#' @param sdf a \code{SigDF}
#' @param pval.threshold minimum p-value to mask
#' @param return.pval whether to return p-values, instead of a
#' masked \code{SigDF}
#' @return a \code{SigDF}, or a p-value vector if return.pval is TRUE
#' @examples
#' sdf <- sesameDataGet("HM450.1.TCGA.PAAD")$sdf
#' sum(mask(sdf))
#' sdf_with_mask <- detectionPoobEcdf(sdf)
#' sum(mask(sdf_with_mask))
#' @export
detectionPoobEcdf2 <- function(sdf, return.pval = FALSE, pval.threshold=0.05){

    stopifnot(is(sdf, "SigDF"))

    func <- ecdf(c(
        with(IR(sdf), c(MG,UG)),
        with(IG(sdf), c(MR,UR))))

    ## p-value is the minimium detection p-value of the 2 alleles
    pvals = setNames(pmin(
        with(sdf, 1-func(pmax(MR, UR, na.rm=TRUE))),
        with(sdf, 1-func(pmax(MG, UG, na.rm=TRUE)))), sdf$Probe_ID)
        
    if (return.pval) { return(pvals) }

    addMask(sdf, pvals > pval.threshold)
}

#' @rdname detectionPoobEcdf
#' @export
#' @examples
#' sdf <- sesameDataGet("HM450.1.TCGA.PAAD")$sdf
#' sum(mask(sdf))
#' sdf_with_mask <- pOOBAH(sdf)
#' sum(mask(sdf_with_mask))
pOOBAH <- detectionPoobEcdf

#' @rdname detectionPoobEcdf2
#' @export
#' @examples
#' sdf <- sesameDataGet("HM450.1.TCGA.PAAD")$sdf
#' sum(mask(sdf))
#' sdf_with_mask <- pOOBAH2(sdf)
#' sum(mask(sdf_with_mask))
pOOBAH2 <- detectionPoobEcdf2
