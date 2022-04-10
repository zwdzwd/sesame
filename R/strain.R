mouseBetaToAF <- function(betas) {
    se <- sesameDataGet('MM285.addressStrain')$strain_snps
    rd <- rowData(se)
    af <- betas[rd$Probe_ID]
    af[rd$flipToAF] <- 1 - af[rd$flipToAF]
    af
}

#' Infer strain information for mouse array
#'
#' @param sdf SigDF
#' @param min_frac_dt minimum fraction of detected signal (DEFAULT: 0.2)
#' otherwise, we give up strain inference and return NA.
#' @param return.probability return probability vector for all strains
#' @param return.pval return p-value
#' @param return.strain return strain name
#' @param verbose print more messages
#' @return a list of best guess, p-value of the best guess
#' and the probabilities of all strains
#' @examples
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet('MM285.1.SigDF')
#' inferStrain(sdf, return.strain = TRUE)
#' sdf.strain <- inferStrain(sdf)
#' @import tibble
#' @export
inferStrain <- function(
    sdf, return.strain = FALSE, return.probability = FALSE,
    return.pval = FALSE, min_frac_dt = 0.2, verbose = FALSE) {

    addr <- sesameDataGet("MM285.addressStrain")
    se <- addr$strain_snps
    cd <- SummarizedExperiment::colData(se)
    rd <- SummarizedExperiment::rowData(se)
    md <- metadata(se)
    
    ## C57BL_6J is the first strain in the table
    strain_snps <- rd[,which(colnames(rd)=="C57BL_6J"):ncol(rd)]

    ## give up if the success rate is low
    pvals <- pOOBAH(sdf, return.pval=TRUE)
    if (sum(pvals[rd$Probe_ID] < 0.05) / nrow(rd) < min_frac_dt) {
        if (return.strain) { return(NA)
        } else if (return.probability) { return(rep(NA, ncol(strain_snps)))
        } else if (return.pval) { return(NA)
        } else { return(sdfMsg(sdf, verbose,
            "Abort strain inference for low detection rate.")) }
    }

    vafs <- getBetas(dyeBiasNL(noob(sdf)), mask=FALSE)[rd$Probe_ID]
    vafs[is.na(vafs)] <- 0.5 # just in case
    vafs[rd$flipToAF] <- 1 - vafs[rd$flipToAF]
    
    probes <- intersect(names(vafs), rd$Probe_ID[rd$QC!="FAIL"])
    vafs <- vafs[probes]
    bbloglik <- vapply(strain_snps[match(probes, rd$Probe_ID),],
        function(x) sum(log(dnorm(x - vafs, mean=0, sd=0.8))), numeric(1))
    probs <- setNames(exp(bbloglik - max(bbloglik)), colnames(strain_snps))

    best.index <- which.max(probs)
    strain <- names(best.index)
    if (return.strain) {
        strain # addr$strain[[strain]][c("JAX_ID","MGP_ID")]
    } else if (return.probability) {
        probs / sum(probs)
    } else if (return.pval) {
        1 - probs[best.index] / sum(probs)
    } else {
        updateSigDF(sdf, strain = strain, addr = addr, verbose = verbose) }
}

#' Compare Strain SNPs with a reference panel
#'
#' @param betas beta value vector or matrix (for multiple samples)
#' @param show_sample_names whether to show sample name
#' @param query_width optional argument for adjusting query width
#' @return grid object that contrast the target sample with
#' pre-built mouse strain reference
#' @importFrom S4Vectors metadata
#' @import wheatmap
#' @export
#' @examples
#' sesameDataCache() # if not done yet
#' compareMouseStrainReference()
#' @export
compareMouseStrainReference <- function(
    betas = NULL, show_sample_names = FALSE, query_width = NULL) {

    ## betas = NULL; show_sample_names = FALSE;
    se <- sesameDataGet("MM285.addressStrain")$strain_snps
    
    cd <- as_tibble(SummarizedExperiment::colData(se))
    rd <- as_tibble(SummarizedExperiment::rowData(se))
    md <- metadata(se)
    se <- se[rd$QC != "FAIL",]; rd <- rd[rd$QC != "FAIL",]

    if (!is.null(betas) && is.null(dim(betas))) { # in case a vector
        betas <- cbind(betas) }

    afs <- do.call(rbind, lapply(seq_along(rd$flipToAF), function(i)
        if(xor(rd$flipToAF[i], rd$flipForRefBias[i])) {
            1-assay(se)[i,]} else {assay(se)[i,]}))
    rownames(afs) <- rd$Probe_ID

    stops <- c("white", "black")
    g <- WHeatmap(afs, cmp=CMPar(stop.points=stops, dmin=0, dmax=1),
        xticklabels = show_sample_names, xticklabels.n=ncol(afs), name="b1")
    if (!is.null(betas)) {          # query samples
        afs2 <- do.call(rbind, lapply(seq_along(rd$flipToAF), function(i) {
            if(xor(rd$flipToAF[i], rd$flipForRefBias[i])) {
                1 - betas[rd$Probe_ID[i],]
            } else { betas[rd$Probe_ID[i],] }}))
        g <- g + WHeatmap(afs2, RightOf("b1", width=query_width),
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
        height=0.03)
    ## g <- g + WLegendV('bh', Beneath(pad=0.06))
    g + WCustomize(mar.bottom=0.15, mar.right=0.06)
}
