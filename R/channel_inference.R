#' Infer and reset color channel for Type-I probes instead of
#' using what is specified in manifest. The results are stored to
#' sdf@extra$IGG and sdf@extra$IRR slot.
#'
#' IGG => Type-I green that is inferred to be green
#' IRR => Type-I red that is inferred to be red 
#' 
#' @param sdf a \code{SigDF}
#' @param verbose whether to print correction summary
#' @param switch_failed whether to switch failed probes (default to FALSE)
#' @param summary return summarized numbers only.
#' @importFrom matrixStats rowMaxs
#' @importFrom matrixStats rowMins
#' @return a \code{SigDF}, or numerics if summary == TRUE
#' @examples
#'
#' sdf <- sesameDataGet('EPIC.1.LNCaP')$sdf
#' inferTypeIChannel(sdf)
#' 
#' @export
inferTypeIChannel <- function(
    sdf, switch_failed = FALSE, verbose = FALSE, summary = FALSE) {
    
    red_channel <- rbind(IR(sdf), oobR(sdf))
    grn_channel <- rbind(oobG(sdf), IG(sdf))
    n_red <- nrow(IR(sdf))
    red_idx0 <- seq_len(nrow(red_channel)) <= nrow(IR(sdf)) # old red index

    ## If there are NA in the probe intensity, exclude these probes.
    ## This is rare and usually occurred when manifest is not complete
    no_na <- complete.cases(cbind(red_channel, grn_channel))
    if (!all(no_na)) {
        red_channel <- red_channel[no_na,]
        grn_channel <- grn_channel[no_na,]
        red_idx0 <- red_idx0[no_na]
        if (verbose) {
            message(
                'Warning! ', sum(!no_na),
                ' Infinium I probes are excluded for having NA intensity.')
        }
    }
    
    red_max <- rowMaxs(red_channel)
    grn_max <- rowMaxs(grn_channel)
    red_idx <- red_max > grn_max # new red index


    inf1_idx = which(sdf$col != "2")
    sdf1 = sdf[inf1_idx,]
    red_max = with(sdf1, pmax(MR, UR))
    grn_max = with(sdf1, pmax(MG, UG))
    new_col = ifelse(red_max > grn_max, "R", "G")
    bg_max = quantile(c(
        with(sdf1[new_col == "R"], c(MG,UG)),
        with(sdf1[new_col == "G"], c(MR,UR))), 0.95)
    
    if (!switch_failed) {
        idx = pmax(red_max, grn_max) < bg_max
        new_col[idx] = sdf1$col[idx]
    }
    sdf$col[inf1_idx] = new_col
    

    ## ## stop inference when in-band signal is lower than a minimum
    ## min_ib <- quantile(
    ##     pmin(rowMins(red_channel), rowMins(grn_channel)), 0.95)
    
    ## big_idx <- pmax(red_max, grn_max) > min_ib # in-band is big enough?

    smry <- c(
        R2R = sum(sdf1$col == "R" & new_col == "R"),
        G2G = sum(sdf1$col == "G" & new_col == "G"),
        R2G = sum(sdf1$col == "R" & new_col == "G"),
        G2R = sum(sdf1$col == "G" & new_col == "R"))
    
    ## if (!switch_failed)
    ##     red_idx <- ifelse(big_idx, red_idx, red_idx0)
    
    ## sdf@extra[['IRR']] <- red_idx[seq_len(n_red)]
    ## sdf@extra[['IGG']] <- !red_idx[(n_red+1):length(red_idx)]
    
    if (summary) { return(smry) }

    if (verbose) {
        message(
            'Infinium-I color channel reset:\n',
            'R>R: ', smry['R2R'], '\n',
            'G>G: ', smry['G2G'], '\n',
            'R>G: ', smry['R2G'], '\n',
            'G>R: ', smry['G2R'], '\n')
    }

    sdf
}

