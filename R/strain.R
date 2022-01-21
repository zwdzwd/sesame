## mouseBetaToAF <- function(betas) {

##     mft <- sesameDataGet('MM285.mm10.manifest')
##     mft <- mft[grep('^rs', names(mft))]
##     ## flip AF based on manifest annotation
##     design <- GenomicRanges::mcols(mft)[['design']]
##     toFlip <- !setNames(as.logical(substr(
##         design, nchar(design), nchar(design))), names(mft))

##     vafs <- betas[grep('^rs', names(betas))]
##     vafs[toFlip[names(vafs)]] <- 1-vafs[toFlip[names(vafs)]]
##     vafs
## }

#' Infer strain information for mouse array
#'
#' @param sdf SigDF
#' @param addr addressStrain, infer and download if not given
#' @return a list of best guess, p-value of the best guess
#' and the probabilities of all strains
#' @examples
#' sesameDataCache("MM285") # if not done yet
#' sdf <- sesameDataGet('MM285.1.SigDF')
#' inferStrain(sdf, return.strain = TRUE)
#' sdf.strain <- inferStrain(sdf)
#' @import tibble
#' @export
inferStrain <- function(
    sdf, addr = NULL,
    return.strain = FALSE, return.probability = FALSE, return.pval = FALSE) {

    ## sdf <- sesameDataGet('MM285.1.SigDF'); betas <- getBetas(dyeBiasNL(noob(sdf)))
    ## vafs <- mouseBetaToAF(betas)
    ## vafs[is.na(vafs)] <- 0.5 # impute vaf if missing
    ## if (is.null(strain_snp_table)) {
    ##     ## TODO: use MM285.strainSNPs
    ##     strain_snp_table <- sesameDataGet('MM285.strain.snp.table')
    ## }
    
    addr <- sesameDataGet("MM285.addressStrain")
    se <- addr$strain_snps
    cd <- SummarizedExperiment::colData(se)
    rd <- SummarizedExperiment::rowData(se)
    md <- metadata(se)
    
    ## C57BL_6J is the first strain in the table
    strain_snps <- rd[,which(colnames(rd)=="C57BL_6J"):ncol(rd)]

    betas <- getBetas(dyeBiasNL(noob(sdf)), mask=FALSE)
    vafs <- betas[rd$Probe_ID]
    vafs[is.na(vafs)] <- 0.5 # just in case
    vafs[rd$flipToAF] <- 1 - vafs[rd$flipToAF]
    
    probes <- intersect(names(vafs), rd$Probe_ID[rd$QC!="FAIL"])
    vafs <- vafs[probes]
    bbloglik <- vapply(strain_snps[match(probes, rd$Probe_ID),],
        function(x) sum(log(dnorm(x - vafs, mean=0, sd=0.8))), numeric(1))
    ## bb <- vapply(probes, function(p) {
    ##     ## vafs[p]; head(sort(setNames(dnorm(vafs[p], mean=as.numeric(strain_snps[match(p, rd$Probe_ID),]), sd=0.8), colnames(strain_snps)), decreasing=T),n=5)
    ##     dnorm(vafs[p], mean=as.numeric(strain_snps[match(p, rd$Probe_ID),]), sd=0.8)
    ## }, numeric(ncol(strain_snps)))
    ## bbloglik <- apply(bb,1,function(x) sum(log(x),na.rm=TRUE))
    probs <- setNames(exp(bbloglik - max(bbloglik)), colnames(strain_snps))
    best.index <- which.max(probs)

    strain <- names(best.index)
    if (return.strain) {
        strain
    } else if (return.probability) {
        probs / sum(probs)
    } else if (return.pval) {
        1 - probs[best.index] / sum(probs)
    } else {
        updateSigDF(sdf, strain = strain, addr = addr)
    }
}

#' Compare Strain SNPs with a reference panel
#'
#' @param betas beta value vector or matrix (for multiple samples)
#' @param show_sample_names whether to show sample name
#' @return grid object that contrast the target sample with
#' pre-built mouse strain reference
#' @import wheatmap
#' @importFrom S4Vectors metadata
#' @export
#' @examples
#' sesameDataCache("MM285") # if not done yet
#' compareMouseStrainReference()
#' @export
compareMouseStrainReference <- function(
    betas = NULL, show_sample_names = FALSE) {

    ## betas = NULL; show_sample_names = FALSE;
    se = sesameDataGet("MM285.addressStrain")$strain_snps
    pkgTest("wheatmap")

    cd <- as_tibble(SummarizedExperiment::colData(se))
    rd <- as_tibble(SummarizedExperiment::rowData(se))
    md <- metadata(se)
    if (!is.null(betas) && is.null(dim(betas))) { # in case a vector
        betas <- cbind(betas)
    }

    se <- se[rd$QC != "FAIL",]
    rd <- rd[rd$QC != "FAIL",]

    afs <- do.call(rbind, lapply(seq_along(rd$flipToAF), function(i)
        if(xor(rd$flipToAF[i], rd$flipForRefBias[i])) {
            1-assay(se)[i,]} else {assay(se)[i,]}))
    rownames(afs) <- rd$Probe_ID

    stops <- c("white", "black")
    g <- WHeatmap(afs, cmp=CMPar(stop.points=stops, dmin=0, dmax=1),
        xticklabels = show_sample_names, xticklabels.n=ncol(afs), name="b1")
    if (!is.null(betas)) {          # query samples
        afs2 <- do.call(rbind, lapply(
            seq_along(rd$flipToAF), function(i) {
                if(xor(rd$flipToAF[i], rd$flipForRefBias[i])) {
                    1 - betas[rd$Probe_ID[i],]
                } else {
                    betas[rd$Probe_ID[i],]
                }}))
        g <- g + WHeatmap(afs2, RightOf("b1"),
            cmp=CMPar(stop.points=stops, dmin=0, dmax=1),
            name="b2", xticklabels=TRUE, xticklabels.n=ncol(betas))
        right <- "b2"
    } else { # in case target is not given, plot just the reference
        right <- "b1"
    }
    
    ## branch color bar (vertical)
    g <- g + WColorBarV(rd$BranchLong, RightOf(right, width=0.03),
        cmp=CMPar(label2color=md$strain.colors), name="bh")
    ## strain color bar (horizontal)
    g <- g + WColorBarH(cd$strain, TopOf("b1",height=0.03),
        cmp=CMPar(label2color=md$strain.colors), name="st")
    ## legends
    g <- g + WLegendV("st",
        TopRightOf("bh", just=c('left','top'), h.pad=0.02),
        height=0.02)
    ## g <- g + WLegendV('bh', Beneath(pad=0.06))
    g + WCustomize(mar.bottom=0.15, mar.right=0.06)
}
