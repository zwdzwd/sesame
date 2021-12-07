mouseBetaToAF = function(betas) {

    mft = sesameDataGet('MM285.mm10.manifest')
    mft = mft[grep('^rs', names(mft))]
    ## flip AF based on manifest annotation
    design = GenomicRanges::mcols(mft)[['design']]
    toFlip = !setNames(as.logical(substr(
        design, nchar(design), nchar(design))), names(mft))

    vafs = betas[grep('^rs', names(betas))]
    vafs[toFlip[names(vafs)]] = 1-vafs[toFlip[names(vafs)]]
    vafs
}

#' Infer strain information for mouse array
#'
#' @param betas beta value vector from which VAFs are extracted
#' @param strain_snp_table if not given download the default from sesameData
#' @return a list of best guess, p-value of the best guess
#' and the probabilities of all strains
#' @examples
#' sesameDataCache("MM285") # if not done yet
#' sdf = sesameDataGet('MM285.1.SigDF')
#' betas = getBetas(dyeBiasNL(noob(sdf)))
#' inferStrain(betas)
#' @import tibble
#' @export
inferStrain <- function(betas, strain_snp_table = NULL) {

    vafs = mouseBetaToAF(betas)

    if (is.null(strain_snp_table)) {
        ## TODO: use MM285.strainSNPs
        strain_snp_table <- sesameDataGet('MM285.strain.snp.table')
    }
    vafs[is.na(vafs)] <- 0.5
    probes <- intersect(names(vafs), rownames(strain_snp_table))
    bb <- vapply(probes, function(p) {
        dnorm(vafs[p], mean=strain_snp_table[p,], sd=0.8)
    }, numeric(ncol(strain_snp_table)))
    bbloglik <- apply(bb,1,function(x) sum(log(x),na.rm=TRUE))
    probs <- exp(bbloglik - max(bbloglik))
    best.index <- which.max(probs)

    list(
        best = names(best.index),
        pval = (sum(probs) - probs[best.index]) / sum(probs),
        probs = probs/sum(probs))
}

#' Compare Strain SNPs with a reference panel
#'
#' @param betas beta value vector or matrix (for multiple samples)
#' @param show_sample_names whether to show sample name
#' @return grid object that contrast the target sample with
#' pre-built mouse strain reference
#' @import wheatmap
#' @import S4Vectors
#' @export
#' @examples
#' sesameDataCache("MM285") # if not done yet
#' compareMouseStrainReference()
#' @export
compareMouseStrainReference = function(
    betas = NULL, show_sample_names = FALSE) {

    se = sesameDataGet("MM285.strainSNPs") # TODO
    pkgTest("wheatmap")

    cd = as_tibble(SummarizedExperiment::colData(se))
    rd = as_tibble(SummarizedExperiment::rowData(se))
    md = metadata(se)
    if (!is.null(betas) && is.null(dim(betas))) { # in case a vector
        betas = cbind(betas)
    }

    afs = do.call(rbind, lapply(seq_along(rd$flipForRefBias), function(i)
        if(rd$flipForRefBias[i]) {1-assay(se)[i,]} else {assay(se)[i,]}))
    rownames(afs) = rd$Probe_ID

    stops = c("white", "black")
    g = WHeatmap(afs, cmp=CMPar(stop.points=stops, dmin=0, dmax=1),
        xticklabels = show_sample_names, xticklabels.n=ncol(afs), name="b1")
    if (!is.null(betas)) {          # query samples
        afs2 = do.call(rbind, lapply(seq_along(rd$flipForRefBias), function(i) {
            if(rd$flipForRefBias[i]) {
                1 - betas[rd$Probe_ID[i],]
            } else {
                betas[rd$Probe_ID[i],]
            }}))
        g = g + WHeatmap(afs2, RightOf("b1"),
            cmp=CMPar(stop.points=stops, dmin=0, dmax=1),
            name="b2", xticklabels=TRUE, xticklabels.n=ncol(betas))
        right = "b2"
    } else { # in case target is not given, plot just the reference
        right = "b1"
    }
    
    ## branch color bar (vertical)
    g = g + WColorBarV(rd$Branch, RightOf(right, width=0.03),
        cmp=CMPar(label2color=md$branch_color), name="bh")
    ## strain color bar (horizontal)
    g = g + WColorBarH(cd$strain, TopOf("b1",height=0.03),
        cmp=CMPar(label2color=md$strain_color), name="st")
    ## legends
    g = g + WLegendV("st", TopRightOf("bh", just=c('left','top'), h.pad=0.02),
        height=0.02)
    g = g + WLegendV('bh', Beneath(pad=0.06))
    g + WCustomize(mar.bottom=0.15, mar.right=0.06)
}
