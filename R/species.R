speciesInfo <- function(addr, species) {
    res <- addr$species[[species]]
    res[c("scientificName", "taxonID", "commonName", "assembly")]
}

#' Set color and mask using strain/species-specific manifest
#'
#' also sets attr(,"species")
#'
#' @param sdf a \code{SigDF}
#' @param species the species the sample is considered to be
#' @param strain the strain the sample is considered to be
#' @param addr species-specific address species, optional
#' @return a \code{SigDF} with updated color channel and mask
#' @examples
#' sdf <- sesameDataGet('Mammal40.1.SigDF')
#' sdf_mouse <- updateSigDF(sdf, species="mus_musculus")
#' 
#' @export
updateSigDF <- function(sdf, species = NULL, strain = NULL, addr = NULL) {

    if (!is.null(species)) {
        if (is.null(addr)) {
            addr <- sesameDataGet(sprintf(
                "%s.addressSpecies", sdfPlatform(sdf)))
        }
        stopifnot(species %in% names(addr$species))
        addrS <- addr$species[[species]]
        message(sprintf("Update using species: %s", species))
    } else if (!is.null(strain)) {
        if (is.null(addr)) {
            addr <- sesameDataGet(sprintf("%s.addressStrain", sdfPlatform(sdf)))
        }
        stopifnot(strain %in% names(addr$strain))
        addrS <- addr$strain[[strain]]
        message(sprintf("Update using strain: %s", strain))
    } else {
        stop("Please specify a species or strain.")
    }
    
    ## set color
    m <- match(sdf$Probe_ID, addr$ordering$Probe_ID)
    ## matched Inf-I probes with non-NA value
    ## (NA can be mapping issues)
    m_idx <- (!is.na(m)) & !is.na(addrS$col[m]) & (sdf$col != "2")
    nc <- as.character(addrS$col[m[m_idx]])
    nc[is.na(nc)] <- '2'
    sdf$col[m_idx] <- factor(nc, levels=c("G","R","2"))

    ## set mask, NA is masked by default
    m_idx <- (!is.na(m))
    sdf$mask <- TRUE
    nm <- addrS$mask[m[m_idx]]
    sdf$mask[!nm] <- FALSE
    attr(sdf, "species") <- species
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
#' increase this number can sometimes improve accuracy (DEFAULT: 1000)
#' @param threshold.pos pvalue < threshold.pos are considered positive
#' (default: 0.01).
#' @param threshold.neg pvalue > threshold.neg are considered negative
#' (default: 0.2).
#' @param return.auc return AUC calculated, override return.species
#' @param return.species return a string to represent species
#' @return a SigDF
#' @examples 
#' sdf <- sesameDataGet("MM285.1.SigDF")
#' sdf <- inferSpecies(sdf)
#'
#' ## all available species
#' all_species <- names(sesameDataGet(sprintf(
#'   "%s.addressSpecies", sdfPlatform(sdf)))$species)
#' 
#' @export
inferSpecies <- function(sdf, topN = 1000,
    threshold.pos = 0.01, threshold.neg = 0.1,
    return.auc = FALSE, return.species = FALSE) {

    addr <- sesameDataGet(sprintf("%s.addressSpecies", sdfPlatform(sdf)))
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
    topN1 <- min(length(neg_probes),length(pos_probes), topN)
    pos_probes <- pos_probes[seq_len(topN1)]
    neg_probes <- neg_probes[seq_len(topN1)]
    
    ## for positive probes (pvalue <= 0.01), y_true = 1
    ## for negative probes (pvalue > 0.1), y_true = 0
    y_true <- structure(c(
        rep(TRUE,length(pos_probes)),rep(FALSE,length(neg_probes))),
        names = c(names(pos_probes), names(neg_probes)))
    
    ## No useful signal, use reference
    if (length(y_true) == 0){
        warning("Lack of useful signal. Use reference.")
        species <- addr$reference
        if (return.auc){ return(NULL);
        } else if (return.species) { return(speciesInfo(addr, species));
        } else { return(updateSigDF(sdf, species=species, addr=addr)); }}
    
    ## calculate AUC based on y_true and y_pred
    n1 <- as.numeric(sum(y_true))
    n2 <- as.numeric(sum(!y_true))
    df_as <- df_as[names(y_true),,drop = FALSE]
    ## df_as[df_as < 35] <- 35 # all under 35 is qualitatively the same
    auc <- vapply(colnames(df_as),function(s) {
        R1 <- sum(rank(df_as[,s])[seq_along(pos_probes)])
        U1 <- R1 - n1 * (n1 + 1)/2
        U1/(n1 * n2)}, numeric(1))

    ## the following is a empirical ladder where one is going to call
    ## reference for lack of negative probes
    if (success.rate >= 0.95 || 
            (success.rate >= 0.80 && max(auc) < 0.50)) {            
        message("Lack of negative probes. Use reference.")
        species <- addr$reference
        if (return.auc){ return(auc);
        } else if (return.species) { return(speciesInfo(addr, species));
        } else { return(updateSigDF(sdf, species=species, addr=addr)); }}

    species <- names(which.max(auc))
    if (return.auc) {
        auc
    } else if (return.species) {
        speciesInfo(addr, species)
    } else {
        updateSigDF(sdf, species=species, addr=addr)}
}

#' Map the SDF (from overlap array platforms)
#' Replicates are merged by picking the best detection
#'
#' @param sdf a \code{SigDF} object
#' @return a named numeric vector for beta values
#' @examples
#' sdf <- sesameDataGet("Mammal40.1.SigDF")
#' betas <- mapToMammal40(sdf[1:10,])
#' @export
mapToMammal40 <- function(sdf) {
    addr <- sesameDataGet("Mammal40.address")
    betas <- getBetas(sdf, collapseToPfx = TRUE)[addr$ordering$Probe_ID]
    names(betas) <- addr$ordering$Probe_ID
    betas
}
