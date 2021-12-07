
#' Mask detection by intermediate beta values
#'
#' @param sdf a \code{SigDF}
#' @param return.pval whether to return p-values, instead of a SigDF
#' @param pval.threshold minimum p-value to mask
#' @param capMU the maximum M+U to search for intermediate betas
#' @param window window size for smoothing and beta fraction calc.
#' @return a \code{SigDF} with mask added
#' @examples
#' sdf <- sesameDataGet("EPIC.1.SigDF")
#' sum(sdf$mask)
#' sum(detectionIB(sdf)$mask)
#' @export
detectionIB = function(
    sdf, return.pval = FALSE, pval.threshold = 0.05,
    capMU = 3000, window = 100) {

    df = signalMU(sdf)
    df$MU = df$M + df$U
    df$beta = df$M / (df$M + df$U)
    df = df[order(df$MU),]

    ## fraction of the intermediate beta
    fmid = vapply(seq_len(nrow(df)-window), function(i) {
        sum(abs(df$beta[i:(i+window-1)] - 0.5)<0.1) / window;
    }, numeric(1))

    ## the range to search for intermediate betas
    error_max_index = min(
        which(fmid <= sort(fmid)[1] + 0.01)[1],
        which(df$MU > capMU)[1])
    df_err = df[seq_len(error_max_index),]
    errors = df_err$MU[abs(df_err$beta - 0.5) < 0.2]
    
    pvals = setNames(1-ecdf(errors)(df$MU), df$Probe_ID)
    pvals[is.na(pvals)] = 1.0 # set NA to 1

    if (return.pval) { return(pvals) }
    addMask(sdf, pvals > pval.threshold)
}

negControls <- function(sdf) {
    stopifnot(is(sdf, "SigDF"))

    ## controls from attributes
    negctls = controls(sdf)[grep(
        'negative', tolower(controls(sdf)$type)),]
    if (!is.null(negctls) && nrow(negctls) > 0) {
        stopifnot(all(c("G","R","col") %in% colnames(negctls)))
        negctls = negctls[negctls$col!=-99, c("G","R")]
    } else { # controls from normal probes
        negctls = as.data.frame(
            sdf[grep("ctl-Negative", sdf$Probe_ID), c("UG", "UR")])
        colnames(negctls) = c("G","R")
    }
    
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
#' sdf <- sesameDataGet("EPIC.1.SigDF")
#' sum(sdf$mask)
#' sum(detectionPnegEcdf(sdf)$mask)
#' @import methods
#' @export
detectionPnegEcdf <- function(sdf, return.pval = FALSE, pval.threshold=0.05) {

    stopifnot(is(sdf, "SigDF"))

    negctls <- negControls(sdf)
    funcG <- ecdf(negctls$G)
    funcR <- ecdf(negctls$R)

    ## p-value is the minimium detection p-value of the 2 alleles
    pvals = setNames(pmin(
        1-funcR(pmax(sdf$MR, sdf$UR, na.rm=TRUE)),
        1-funcG(pmax(sdf$MG, sdf$UG, na.rm=TRUE))), sdf$Probe_ID)
        
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
#' @param combine.neg whether to combine negative control probes with
#' the out-of-band probes in simulating the signal background
#' @return a \code{SigDF}, or a p-value vector if return.pval is TRUE
#' @examples
#' sdf <- sesameDataGet("EPIC.1.SigDF")
#' sum(sdf$mask)
#' sum(detectionPoobEcdf(sdf)$mask)
#' 
#' @export
detectionPoobEcdf <- function(sdf, return.pval = FALSE,
    combine.neg = TRUE, pval.threshold=0.05) {

    stopifnot(is(sdf, "SigDF"))

    dG = InfIG(sdf); dR = InfIR(sdf)
    bgG = c(dR$MG, dR$UG)
    bgR = c(dG$MR, dG$UR)
    if (combine.neg) {
        neg = negControls(sdf)
        bgG = c(bgG, neg$G)
        bgR = c(bgR, neg$R)
    }
    funcG = ecdf(bgG)
    funcR = ecdf(bgR)

    ## p-value is the minimium detection p-value of the 2 alleles
    pvals = setNames(pmin(
        1-funcR(pmax(sdf$MR, sdf$UR, na.rm=TRUE)),
        1-funcG(pmax(sdf$MG, sdf$UG, na.rm=TRUE))), sdf$Probe_ID)
        
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
#' @param combine.neg whether to combine negative control probes with
#' the out-of-band probes in simulating the signal background
#' @return a \code{SigDF}, or a p-value vector if return.pval is TRUE
#' @examples
#' sdf <- sesameDataGet("EPIC.1.SigDF")
#' sum(sdf$mask)
#' sdf <- detectionPoobEcdf2(sdf)
#' sum(sdf$mask)
#' @export
detectionPoobEcdf2 <- function(sdf, return.pval = FALSE,
    combine.neg = TRUE, pval.threshold=0.05){

    stopifnot(is(sdf, "SigDF"))

    dG = InfIG(sdf); dR = InfIR(sdf)
    bg = c(dR$MG,dR$UG,dG$MR,dG$UR)
    if (combine.neg) {
        neg = negControls(sdf)
        bg = c(bg, c(neg$G, neg$R))
    }
    func <- ecdf(bg)
    
    ## p-value is the minimium detection p-value of the 2 alleles
    pvals = setNames(pmin(
        1-func(pmax(sdf$MR, sdf$UR, na.rm=TRUE)),
        1-func(pmax(sdf$MG, sdf$UG, na.rm=TRUE))), sdf$Probe_ID)
        
    if (return.pval) { return(pvals) }

    addMask(sdf, pvals > pval.threshold)
}

#' @rdname detectionPoobEcdf
#' @export
#' @examples
#' sdf <- sesameDataGet("EPIC.1.SigDF")
#' sum(sdf$mask)
#' sum(resetMask(sdf)$mask)
#' sum(pOOBAH(sdf, pval.threshold=0.2)$mask)
pOOBAH <- detectionPoobEcdf

#' @rdname detectionPoobEcdf2
#' @export
#' @examples
#' sdf <- sesameDataGet("EPIC.1.SigDF")
#' sum(sdf$mask)
#' sdf <- pOOBAH2(sdf)
#' sum(sdf$mask)
pOOBAH2 <- detectionPoobEcdf2
