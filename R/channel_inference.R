#' Infer and reset color channel for Type-I probes instead of
#' using what is specified in manifest. The results are stored to
#' sset@extra$IGG and sset@extra$IRR slot.
#'
#' IGG => Type-I green that is inferred to be green
#' IRR => Type-I red that is inferred to be red 
#' 
#' @param sset a \code{SigSet}
#' @param verbose whether to print correction summary
#' @param switch_failed whether to switch failed probes (default to FALSE)
#' @param summary return summarized numbers only.
#' @importFrom matrixStats rowMaxs
#' @importFrom matrixStats rowMins
#' @return a \code{SigSet}, or numerics if summary == TRUE
#' @examples
#'
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#' inferTypeIChannel(sset)
#' 
#' @export
inferTypeIChannel <- function(
    sset, switch_failed = FALSE, verbose = FALSE, summary = FALSE) {
    
    red_channel <- rbind(IR(sset), oobR(sset))
    grn_channel <- rbind(oobG(sset), IG(sset))
    n_red <- nrow(IR(sset))
    red_idx0 <- seq_len(nrow(red_channel)) <= nrow(IR(sset)) # old red index

    ## If there are NA in the probe intensity, exclude these probes.
    ## This is rare and usually occurred when manifest is not complete
    no_na <- complete.cases(cbind(red_channel, grn_channel))
    if (!all(no_na)) {
        red_channel <- red_channel[no_na,]
        grn_channel <- grn_channel[no_na,]
        red_idx0 <- red_idx0[no_na]
        if (verbose) {
            message(
                sum(!no_na),
                ' Infinium I probes are excluded for having NA intensity.')
        }
    }
    
    red_max <- rowMaxs(red_channel)
    grn_max <- rowMaxs(grn_channel)
    red_idx <- red_max > grn_max # new red index

    ## stop inference when in-band signal is lower than a minimum
    min_ib <- quantile(
        pmin(rowMins(red_channel), rowMins(grn_channel)), 0.95)
    
    big_idx <- pmax(red_max, grn_max) > min_ib # in-band is big enough?

    smry <- c(
        R2R = sum(red_idx0 & red_idx & big_idx),
        G2G = sum(!red_idx0 & !red_idx & big_idx),
        R2G = sum(red_idx0 & !red_idx & big_idx),
        G2R = sum(!red_idx0 & red_idx & big_idx),
        FailedR = sum(red_idx0 & !big_idx),
        FailedG = sum(!red_idx0 & !big_idx))

    if (!switch_failed)
        red_idx <- ifelse(big_idx, red_idx, red_idx0)
    
    sset@extra[['IRR']] <- red_idx[seq_len(n_red)]
    sset@extra[['IGG']] <- !red_idx[(n_red+1):length(red_idx)]
    
    if (summary) {
        return(smry)
    }

    if (verbose) {
        message(
            'Type-I color channel reset:\n',
            'R>R: ', smry['R2R'], '\n',
            'G>G: ', smry['G2G'], '\n',
            'R>G: ', smry['R2G'], '\n',
            'G>R: ', smry['G2R'], '\n',
            'Red Failed: ', smry['FailedR'], '\n',
            'Grn Failed: ', smry['FailedG'])
    }

    sset
}

## Type-I Grn after correction
IG2 <- function(sset) {
    if (extraHas(sset,'IGG') && extraHas(sset,'IRR')) {
        rbind(sset@IG[sset@extra$IGG,], sset@oobG[!sset@extra$IRR,])
    } else {
        IG(sset)
    }
}

## Type-I Red after correction
IR2 <- function(sset) {
    if (extraHas(sset,'IGG') && extraHas(sset,'IRR')) {
        rbind(sset@IR[sset@extra$IRR,], sset@oobR[!sset@extra$IGG,])
    } else {
        IR(sset)
    }
}

## OOB Grn after correction
oobG2 <- function(sset) {
    if ('IGG' %in% names(sset@extra) && 'IRR' %in% names(sset@extra)) {
        rbind(sset@oobG[sset@extra$IRR,], sset@IG[!sset@extra$IGG,])
    } else { # backward-compatible
        oobG(sset)
    }
}

## OOB Red after correction
oobR2 <- function(sset) {
    if ('IGG' %in% names(sset@extra) && 'IRR' %in% names(sset@extra)) {
        rbind(sset@oobR[sset@extra$IGG,], sset@IR[!sset@extra$IRR,])
    } else { # backward-compatible
        oobR(sset)
    }
}


