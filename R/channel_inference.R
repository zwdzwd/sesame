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
#' @param mask_failed whether to mask failed probes (default to FALSE)
#' @param summary return summarized numbers only.
#' @return a \code{SigDF}, or numerics if summary == TRUE
#' @examples
#'
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#' inferInfiniumIChannel(sdf)
#' 
#' @export
inferInfiniumIChannel <- function(
    sdf, switch_failed = FALSE, mask_failed = FALSE,
    verbose = FALSE, summary = FALSE) {
    
    inf1_idx <- which(sdf$col != "2")
    sdf1 <- sdf[inf1_idx,]
    red_max <- pmax(sdf1$MR, sdf1$UR)
    grn_max <- pmax(sdf1$MG, sdf1$UG)
    new_col <- factor(ifelse(
        red_max > grn_max, "R", "G"), levels=c("G","R","2"))
    d1R <- sdf1[new_col == "R",]
    d1G <- sdf1[new_col == "G",]
    bg_max <- quantile(c(d1R$MG,d1R$UG,d1G$MR,d1G$UR), 0.95, na.rm=TRUE)

    ## revert to the original for failed probes if so desire
    idx <- (is.na(red_max) | is.na(grn_max) | pmax(red_max, grn_max) < bg_max)
    if (!switch_failed) {
        new_col[idx] <- sdf1$col[idx]
    }
    if (mask_failed) {
        sdf$mask[inf1_idx[idx]] <- TRUE
    }
    sdf$col[inf1_idx] <- factor(new_col, levels=c("G","R","2"))

    smry <- c(
        R2R = sum(sdf1$col == "R" & new_col == "R", na.rm=TRUE),
        G2G = sum(sdf1$col == "G" & new_col == "G", na.rm=TRUE),
        R2G = sum(sdf1$col == "R" & new_col == "G", na.rm=TRUE),
        G2R = sum(sdf1$col == "G" & new_col == "R", na.rm=TRUE))
    
    if (summary) { return(smry) }

    sdfMsg(sdf, verbose, "%s: R>R:%d;G>G:%d;R>G:%d;G>R:%d",
        "Infinium-I color channel reset",
        smry["R2R"], smry["G2G"], smry["R2G"], smry["G2R"])
}

