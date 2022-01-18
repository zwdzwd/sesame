speciesInfo <- function(addr, species) {
    res <- addr$species[[species]]
    res[c("scientificName", "taxonID", "commonName", "assembly")]
}

#' Set color and mask using species-specific manifest
#'
#' @param sdf a \code{SigDF}
#' @param species the species the sample is considered to be
#' @param addr species-specific address species, optional
#' @return a \code{SigDF} with updated color channel and mask
#' @examples
#' sdf <- sesameDataGet('Mammal40.1.SigDF')
#' sdf_mouse <- updateSigBySpecies(sdf, "mus_musculus")
#' 
#' @export
updateSigBySpecies <- function(sdf, species, addr = NULL) {
    if (is.null(addr)) {
        addr <- sesameDataGet(sprintf("%s.addressSpecies", sdfPlatform(sdf)))
    }
    stopifnot(species %in% names(addr$species))
    addr.sp <- addr$species[[species]]
    message(sprintf("Update using species: %s", species))
    
    ## set color
    m <- match(sdf$Probe_ID, addr$ordering$Probe_ID)
    ## matched Inf-I probes with non-NA value
    ## (NA can be mapping issues)
    m_idx <- (!is.na(m)) & !is.na(addr.sp$col[m]) & (sdf$col != "2")
    nc <- as.character(addr.sp$col[m[m_idx]])
    nc[is.na(nc)] <- '2'
    sdf$col[m_idx] <- factor(nc, levels=c("G","R","2"))

    ## set mask, NA is masked by default
    m_idx <- (!is.na(m))
    sdf$mask <- TRUE
    nm <- addr.sp$mask[m[m_idx]]
    sdf$mask[!nm] <- FALSE
    sdf
}

#' Infer Species
#'
#' We infer species based on probes pvalues and alignment score.
#' AUC was calculated for each specie, y_true is 1 or 0 
#' for pval < threshold.pos or pval > threshold.neg, respeceively,
#'
#' @param sdf a \code{SigDF}
#' @param topN Top n positive and negative probes used to infer species.
#' @param threshold.pos pvalue < threshold.pos are considered positive
#' (default: 0.01).
#' @param threshold.neg pvalue > threshold.neg are considered negative
#' (default: 0.2).
#' @param return.auc return AUC calculated, override return.species
#' @param return.species return a string to represent species
#' @param threshold.success.rate threshold of success rate to determine
#' mouse species.
#' @return a new SigDF with updated color channel specification and
#' masking unless return.auc or return.species is set to TRUE
#'
#' @examples 
#' if (FALSE) { ## remove this, testing doesn't allow large file caching
#'   sdf <- sesameDataGet("MM285.1.SigDF")
#'   sdf <- inferSpecies(sdf)
#' }
#' @export
inferSpecies <- function(sdf, topN = 1000,
    threshold.pos = 0.01, threshold.neg = 0.1,
    return.auc = FALSE, return.species = FALSE,
    threshold.success.rate = 0.8) {

    ##sdf = sdf; topN=1000; threshold.pos=0.01; threshold.neg=0.1; ret.max=TRUE; balance=TRUE; threshold.success.rate=0.8
    ## TODO remove use_alternative = TRUE
    addr <- sesameDataGet(sprintf("%s.addressSpecies", sdfPlatform(sdf)), use_alternative = TRUE)
    df_as <- do.call(cbind, lapply(addr$species, function(x) x$AS))
    rownames(df_as) <- addr$ordering$Probe_ID
    
    pvalue <- pOOBAH(sdf, return.pval=TRUE)
    ## shared probes
    pvalue <- pvalue[intersect(names(pvalue),rownames(df_as))]
    ## get positive probes (pvalue <= 0.01) and sort ascendingly
    pos_probes <- sort(pvalue[pvalue <= threshold.pos],decreasing=FALSE)
    ## get negative probes (pvalue >= 0.1) and sort descendingly
    neg_probes <- sort(pvalue[pvalue >= threshold.neg],decreasing=TRUE)
    success.rate <- length(pvalue[pvalue<=0.05]) / length(pvalue)
    
    ## keep the same number of positive and negative probes.
    topN <- min(length(neg_probes),length(pos_probes), topN)
    pos_probes <- pos_probes[seq_len(topN)]
    neg_probes <- neg_probes[seq_len(topN)]
    
    ## for positive probes (pvalue <= 0.01), y_true = 1
    ## for negative probes (pvalue > 0.1), y_true = 0
    y_true <- structure(c(
        rep(1,length(pos_probes)),rep(0,length(neg_probes))),
        names = c(names(pos_probes), names(neg_probes)))
    
    ## y_pred is the alignment score.
    df_as <- df_as[c(names(pos_probes), names(neg_probes)),]
    
    ## No useful signal, use reference
    if (length(y_true) == 0){
        warning("Lack of useful signal. Use reference.")
        species <- addr$reference
        if (return.auc){ return(NULL);
        } else if (return.species) { return(speciesInfo(addr, species));
        } else { return(updateSigBySpecies(sdf, species, addr)); }}
    
    ## calculate AUC based on y_true and y_pred
    auc <- vapply(colnames(df_as),function(s) {
        labels <- as.logical(y_true)
        n1 <- as.numeric(sum(labels))
        n2 <- as.numeric(sum(!labels))
        R1 <- sum(rank(df_as[,s])[labels])
        U1 <- R1 - n1 * (n1 + 1)/2
        U1/(n1 * n2)}, numeric(1))

    ## if success rate is high but max(AUC) is low it can be because the target
    ## species has issue with alignemnt score calibration
    if (success.rate >= threshold.success.rate && max(auc) < 0.6) {
        message("Lack of negative probes. Use reference.")
        species <- addr$reference
        if (return.auc){ return(auc);
        } else if (return.species) { return(speciesInfo(addr, species));
        } else { return(updateSigBySpecies(sdf, species, addr)); }}
    
    species <- names(which.max(auc))
    if (return.auc) {
        auc
    } else if (return.species) {
        speciesInfo(addr, species)
    } else {
        updateSigBySpecies(sdf, species, addr)}
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
