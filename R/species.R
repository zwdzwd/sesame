


#' Infer Species
#'
#' We infer species based on probes pvalues and alignment score.
#' AUC was calculated for each specie, y_true is 1 or 0 
#' for pval < threshold.pos or pval > threshold.neg, respeceively,
#'
#' @param sdf a \code{SigSet}
#' @param df_as a data.frame of alignment score for each probe.
#' @param topN Top n positive and negative probes used to infer species.
#' @param threshold.pos pvalue < threshold.pos are considered positive
#' (default: 0.01).
#' @param threshold.neg pvalue > threshold.neg are considered negative
#' (default: 0.2).
#' @param ret.max whether to return the species with maximal AUC.
#' @param balance whether to balance the postive and negative probe
#' number (default: TRUE).
#' @param threshold.sucess.rate threshold of success rate to determine
#' mouse species.
#' @return a list of auc, pvalue, species (NCBI official species names)
#' and taxid.
#'
#' @examples 
#' if (FALSE) { ## remove this, testing doesn't allow large file caching
#'   sdf <- sesameDataGet("MM285.1.SigDF")
#'   inferSpecies(sdf)
#' }
#' @export
inferSpecies <- function(sdf, df_as = NULL, topN = 3000,
    threshold.pos = 0.01, threshold.neg = 0.1, ret.max = TRUE,
    balance = TRUE, threshold.sucess.rate = 0.8) {

    if (is.null(df_as)) {
        ## Load alignment score (df_as) for candidate species
        df_as <- sesameDataGet(sprintf("%s.alignmentScore", sdfPlatform(sdf)))
    }

    pvalue <- pOOBAH(sdf, return.pval=TRUE)
    ## shared probes
    pvalue <- pvalue[intersect(names(pvalue),rownames(df_as))]
    ## get positive probes (pvalue <= 0.01) and sort ascendingly
    pos_probes <- sort(pvalue[pvalue <= threshold.pos],decreasing=FALSE)
    ## get negative probes (pvalue >= 0.1) and sort descendingly
    neg_probes <- sort(pvalue[pvalue >= threshold.neg],decreasing=TRUE)
    success.rate <- length(pvalue[pvalue<=0.05]) / length(pvalue)
    
    ## if success.rate >0.8, return 'mouse'
    ## TODO: same thing decision for Mammal40?
    if (success.rate >= threshold.sucess.rate && sdfPlatform(sdf) == 'MM285') {
        return(list(auc=1,taxid='10090',species='Mus musculus'))
    }
    
    ## balance means keep the same number of positive and negative probes.
    if (balance) {
        topN <- min(length(neg_probes),length(pos_probes))
    }
    if (length(pos_probes) > topN){
        pos_probes <- pos_probes[seq_len(topN)]
    }
    if (length(neg_probes) > topN){
        neg_probes <- neg_probes[seq_len(topN)]
    }
    
    ## for positive probes (pvalue <= 0.01), y_true = 1
    ## for negative probes (pvalue > 0.1), y_true = 0
    y_true <- structure(c(rep(1,length(pos_probes)),rep(0,length(neg_probes))),
        names = c(names(pos_probes), names(neg_probes)))
    
    ## y_pred is the alignment score.
    df_as <- df_as[c(names(pos_probes), names(neg_probes)),]
    
    if (length(y_true) == 0){
        if (ret.max){
            return(list(auc=NA,species=NA,taxid=NA))
        } else {
            return(as.data.frame(t(vapply(colnames(df_as),function() {
                list(auc=NA,species=NA,taxid=NA)
            }, c(numeric(1), character(1), character(1))))))
        }
    }
    
    ## Calculate AUC based on y_true and y_pred (df_as)
    auc <- vapply(colnames(df_as),function(s) {
        labels <- as.logical(y_true)
        n1 <- sum(labels)
        n2 <- sum(!labels)
        R1 <- sum(rank(df_as[,s])[labels])
        U1 <- R1 - n1 * (n1 + 1)/2
        U1/(n1 * n2)
    }, numeric(1))
    
    ## if ret.max, return only the species with the maximal AUC
    ## a vector of species and AUC would be returned.
    if (ret.max){
        l <- as.list(auc[which.max(auc)])
        l['species'] <- unlist(strsplit(names(l)[1],"\\|"))[1]
        l['taxid'] <- unlist(strsplit(names(l)[1],"\\|"))[2]
        names(l)[1] <- 'auc'
        return(l)
    } else { # return all
        return(auc)
    }
}

#' Map the SDF (from overlap array platforms)
#' Replicates are merged by picking the best detection
#'
#' @param sdf a \code{SigDF} object
#' @return a named numeric vector for beta values
#' @examples
#' sdf <- sesameDataGet("MM285.1.SigDF")
#' betas <- mapToMammal40(sdf)
#' @export
mapToMammal40 <- function(sdf) {
    addr <- sesameDataGet("Mammal40.address")
    betas <- getBetas(sdf, collapseToPfx = TRUE)[addr$ordering$Probe_ID]
    names(betas) <- addr$ordering$Probe_ID
    betas
}
