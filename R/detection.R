
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
detectionIB <- function(
    sdf, return.pval = FALSE, pval.threshold = 0.05,
    capMU = 3000, window = 100) {

    df <- signalMU(sdf)
    df$MU <- df$M + df$U
    df$beta <- df$M / (df$M + df$U)
    df <- df[order(df$MU),]

    ## fraction of the intermediate beta
    fmid <- vapply(seq_len(nrow(df)-window), function(i) {
        sum(abs(df$beta[i:(i+window-1)] - 0.5)<0.1) / window;
    }, numeric(1))

    ## the range to search for intermediate betas
    error_max_index <- min(
        which(fmid <= sort(fmid)[1] + 0.01)[1],
        which(df$MU > capMU)[1])
    df_err <- df[seq_len(error_max_index),]
    errors <- df_err$MU[abs(df_err$beta - 0.5) < 0.2]
    
    pvals <- setNames(1-ecdf(errors)(df$MU), df$Probe_ID)
    pvals[is.na(pvals)] <- 1.0 # set NA to 1

    if (return.pval) { return(pvals) }
    addMask(sdf, pvals > pval.threshold)
}

#' get negative control signal
#'
#' @param sdf a SigDF
#' @return a data frame of negative control signals
negControls <- function(sdf) {
    stopifnot(is(sdf, "SigDF"))

    idx <- grep("negative",tolower(sdf$Probe_ID))
    if (length(idx) > 0) {
        negctls <- sdf[idx,c("UG","UR")]
    } else {
        df <- controls(sdf)
        negctls <- df[grep('negative', tolower(df$Type)), c("UG","UR")]
    }
    colnames(negctls) <- c("G","R")
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
#' @export
detectionPnegEcdf <- function(sdf, return.pval = FALSE, pval.threshold=0.05) {

    stopifnot(is(sdf, "SigDF"))

    negctls <- negControls(sdf)
    funcG <- ecdf(negctls$G)
    funcR <- ecdf(negctls$R)

    ## p-value is the minimium detection p-value of the 2 alleles
    pvals <- setNames(pmin(
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
#' @name pOOBAH
#' @param sdf a \code{SigDF}
#' @param pval.threshold minimum p-value to mask
#' @param return.pval whether to return p-values, instead of a
#' masked \code{SigDF}
#' @param combine.neg whether to combine negative control probes with
#' the out-of-band probes in simulating the signal background
#' @param verbose print more messages
#' @return a \code{SigDF}, or a p-value vector if return.pval is TRUE
#' @examples
#' sdf <- sesameDataGet("EPIC.1.SigDF")
#' sum(sdf$mask)
#' sum(pOOBAH(sdf)$mask)
#' 
#' @export
pOOBAH <- function(sdf, return.pval = FALSE,
    combine.neg = TRUE, pval.threshold=0.05, verbose = FALSE) {

    stopifnot(is(sdf, "SigDF"))
    nmk <- sdf[!(sdf$Probe_ID %in% nonuniqMask(sdfPlatform(sdf))),]
    bgG <- oobG(nmk)
    bgR <- oobR(nmk)
    if (combine.neg) {
        neg <- negControls(sdf)
        bgG <- c(bgG, neg$G)
        bgR <- c(bgR, neg$R)
    }

    ## if there is no background signal, use an empirical prior
    if (sum(!is.na(bgG)) <= 100) { bgG <- seq_len(1000) }
    if (sum(!is.na(bgR)) <= 100) { bgR <- seq_len(1000) }
    
    funcG <- ecdf(bgG)
    funcR <- ecdf(bgR)

    ## p-value is the minimium detection p-value of the 2 alleles
    ## the order is preserved
    pvals <- setNames(pmin(
        1-funcR(pmax(sdf$MR, sdf$UR, na.rm=TRUE)),
        1-funcG(pmax(sdf$MG, sdf$UG, na.rm=TRUE))), sdf$Probe_ID)
        
    if (return.pval) { return(pvals) }

    addMask(sdf, pvals > pval.threshold)
}


## ## if the background is extremely low, be conservative
## bgG95 <- quantile(bgG, 0.95, na.rm=TRUE)
## if (bgG95 < 100) { bgG <- bgG * 100 / bgG95 }
## bgR95 <- quantile(bgR, 0.95, na.rm=TRUE)
## if (bgR95 < 100) { bgR <- bgR * 100 / bgR95 }
