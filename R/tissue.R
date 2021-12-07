reference_plot_se = function(
    betas, se, color=c("blueYellow","jet"), show_sample_names=FALSE) {

    ## top N probes ordered by delta_beta, will
    ## be dominated by certain tissue otherwise
    rd = as.data.frame(rowData(se))
    topN = do.call(c, lapply(split(rd, rd$branch), function(x) {
        x$Probe_ID[order(x$delta_beta)][seq_len(min(nrow(x), 200))] }))
    se = se[rd$Probe_ID %in% topN,]
    
    pkgTest("wheatmap")
    color = match.arg(color)
    if (color == "blueYellow") {
        stop.points = c("blue","yellow")
    } else {
        stop.points = NULL
    }

    cd = as_tibble(colData(se))
    rd = as_tibble(rowData(se))
    md = metadata(se)
    if (!is.null(betas) && is.null(dim(betas))) { # in case a vector
        betas = cbind(betas)
    }

    g = WHeatmap(assay(se), cmp=CMPar(stop.points=stop.points,
        dmin=0, dmax=1), xticklabels = show_sample_names, name="b1") # reference
    if (!is.null(betas)) {          # query samples
        g = g + WHeatmap(betas[rd$Probe_ID,], RightOf("b1"),
            cmp=CMPar(stop.points=stop.points, dmin=0, dmax=1),
            name="b2", xticklabels=TRUE, xticklabels.n=ncol(betas))
        right = "b2"
    } else { # in case target is not given, plot just the reference
        right = "b1"
    }
    ## branch color bar (vertical)
    g = g + WColorBarV(rd$branch, RightOf(right, width=0.03),
        cmp=CMPar(label2color=md$branch_color), name="bh")
    ## tissue color bar (horizontal)
    g = g + WColorBarH(cd$tissue, TopOf("b1",height=0.03),
        cmp=CMPar(label2color=md$tissue_color), name="ti")
    ## legends
    g = g + WLegendV("ti", TopRightOf("bh", just=c('left','top'), h.pad=0.02),
        height=0.02)
    g = g + WLegendV('bh', Beneath(pad=0.06))
    g + WCustomize(mar.bottom=0.15, mar.right=0.06)
}

#' Compare mouse array data with mouse tissue references
#'
#' @param betas matrix of betas for the target sample
#' @param color either blueYellow or fullJet
#' @return grid object that contrast the target sample with
#' pre-built mouse tissue reference
#' @import wheatmap
#' @export
#' @examples
#' 
#' sesameDataCache("MM285") # if not done yet
#' compareMouseTissueReference()
#' sesameDataClearCache()
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
compareMouseTissueReference = function(betas=NULL, color="blueYellow") {
    se = sesameDataGet("MM285.tissueSignature")

    reference_plot_se(betas, se, color=color)
}

#' inferTissue1 infers the tissue of a single sample (as identified through 
#' the branchIDs in the row data of the reference) by reporting independent 
#' composition through cell type deconvolution.
#'
#' @param betas Named vector with probes and their corresponding beta value 
#' measurement
#' @param reference Summarized Experiment with either hypomethylated or 
#' hypermethylated probe selection (row data), sample selection (column data), 
#' meta data, and the betas (assay)
#' @param platform String representing the array type of the betas and 
#' reference
#' @param abs_delta_beta_min Numerical value indicating the absolute minimum 
#' required delta beta for the probe selection criteria
#' @param auc_min Numeric value corresponding to the minimum AUC value 
#' required for a probe to be considered
#' @param coverage_min Numeric value corresponding to the minimum coverage 
#' requirement for a probe to be considered. Coverage is defined here as the 
#' proportion of samples without an NA value at a given probe.
#' @param topN number of probes to at most use for each branch
#'
#' @return inferred tissue as a string
#' @examples
#' sesameDataCache("MM285") # if not done yet
#' sdf = sesameDataGet("MM285.1.SigDF")
#' inferTissue(getBetas(dyeBiasNL(noob(sdf))))
#'
#' sesameDataClearCache()
#'
#' @export
inferTissue = function(betas, reference = NULL, platform = NULL,
    abs_delta_beta_min = 0.3, auc_min = 0.99, coverage_min = 0.80, topN = 15) {

    stopifnot(is.numeric(betas))

    if (is.null(reference)) {
        if (is.null(platform)) {
            platform = inferPlatformFromProbeIDs(names(betas))
        }
        stopifnot(platform %in% c("MM285")) # TODO: add human
        reference = sesameDataGet(sprintf("%s.tissueSignature", platform))
    }
    
    rd = rowData(reference)
    fracs = sort(vapply(unique(rd$branch), function(branch) {
        rd1 = rd[
            rd$branch == branch & abs(rd$delta_beta) >= abs_delta_beta_min, ]

        rd1 = head(rd1[order(-abs(rd1$delta_beta)), ], n = topN)
        
        fracs1 = c(1 - betas[rd1[rd1$delta_beta < 0, "Probe_ID"]],
            betas[rd1[rd1$delta_beta > 0, "Probe_ID"]])

        mean(fracs1, na.rm = TRUE)
    }, numeric(1)), decreasing = TRUE)
    sprintf("[%s](%1.1f) [%s](%1.1f)",
        names(fracs)[1], fracs[1], names(fracs)[2], fracs[2])
    
    ## results = results[!(names(results) %in% ignore_branches)]
    ## cd = meta[match(colnames(results), meta$betas),]
    ## se = SummarizedExperiment(assays=list(results=results), colData=cd)
    ## metadata(se)$tissue_color = metadata(reference)$tissue_color
    ## metadata(se)$branchID_color = metadata(reference)$branchID_color
    ## se
}


