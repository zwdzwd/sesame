#' Compare array data with references (e.g., tissue, cell types)
#'
#' @param ref the reference beta values in SummarizedExperiment.
#' One can download them from the sesameData package. See examples.
#' @param betas matrix of betas for the target sample
#' This argument is optional. If not given, only the reference will be shown.
#' @param stop.points stop points for the color palette.
#' Default to blue, yellow.
#' @param query_width the width of the query beta value matrix
#' @param show_sample_names whether to show sample names (default: FALSE)
#' @return grid object that contrast the target sample with
#' references.
#' @export
#' @examples
#' 
#' sesameDataCache() # if not done yet
#' compareReference(sesameDataGet("MM285.tissueSignature"))
#' sesameDataGet_resetEnv()
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
compareReference <- function(
    ref, betas = NULL, stop.points = NULL, query_width=0.3,
    show_sample_names = FALSE) {

    if (is.null(stop.points)) { stop.points <- c("blue","yellow") }
    
    cd <- as_tibble(colData(ref))
    rd <- as_tibble(rowData(ref))
    md <- metadata(ref)
    if (!is.null(betas) && is.null(dim(betas))) { # in case a vector
        betas <- cbind(betas)
    }

    ## reference
    g <- WHeatmap(assay(ref), cmp=CMPar(stop.points=stop.points,
        dmin=0, dmax=1), xticklabels = show_sample_names, name="b1")
    ## query samples
    if (!is.null(betas)) {
        g <- g + WHeatmap(betas[rd$Probe_ID,], RightOf("b1", width=query_width),
            cmp=CMPar(stop.points=stop.points, dmin=0, dmax=1),
            name="b2", xticklabels = show_sample_names,
            xticklabels.n=ncol(betas))
        right <- "b2"
    } else { # in case target is not given, plot just the reference
        right <- "b1"
    }
    ## branch color bar (vertical)
    g <- g + WColorBarV(rd$branch, RightOf(right, width=0.03),
        cmp=CMPar(label2color=md$branch_color), name="bh")
    ## tissue color bar (horizontal), branch should be replaced by CellType
    g <- g + WColorBarH(cd$branch, TopOf("b1",height=0.03),
        cmp=CMPar(label2color=md$branch_color), name="ti")
    ## legends
    g <- g + WLegendV("ti", TopRightOf("bh", just=c('left','top'), h.pad=0.02),
        height=0.02)
    g + WCustomize(mar.bottom=0.15, mar.right=0.06)
}

#' Compare mouse array data with mouse tissue references
#'
#' @param betas matrix of betas for the target sample
#' This argument is optional. If not given, only the reference will be shown.
#' @param ref the reference beta values in SummarizedExperiment.
#' This argument is optional. If not given, the reference will be downloaded
#' from the sesameData package.
#' @param color either blueYellow or fullJet
#' @param query_width the width of the query beta value matrix
#' @return grid object that contrast the target sample with
#' pre-built mouse tissue reference
#' @export
#' @examples
#' cat("Deprecated, see compareReference")
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
compareMouseTissueReference <- function(
    betas=NULL, ref=NULL, color="blueYellow", query_width=0.3) {
    .Deprecated("compareReference")
}

#' inferTissue infers the tissue of a single sample (as identified through 
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
#' sesameDataCache() # if not done yet
#' sdf <- sesameDataGet("MM285.1.SigDF")
#' inferTissue(getBetas(dyeBiasNL(noob(sdf))))
#'
#' sesameDataGet_resetEnv()
#'
#' @export
inferTissue <- function(betas, reference = NULL, platform = NULL,
    abs_delta_beta_min = 0.3, auc_min = 0.99, coverage_min = 0.80, topN = 15) {

    stopifnot(is.numeric(betas))

    if (is.null(reference)) {
        if (is.null(platform)) {
            platform <- inferPlatformFromProbeIDs(names(betas))
        }
        stopifnot(platform %in% c("MM285")) # TODO: add human
        reference <- sesameDataGet(sprintf("%s.tissueSignature", platform))
    }
    
    rd <- rowData(reference)
    fracs <- sort(vapply(unique(rd$branch), function(branch) {
        rd1 <- rd[
            rd$branch == branch & abs(rd$delta_beta) >= abs_delta_beta_min, ]

        rd1 <- head(rd1[order(-abs(rd1$delta_beta)), ], n = topN)
        
        fracs1 <- c(1 - betas[rd1[rd1$delta_beta < 0, "Probe_ID"]],
            betas[rd1[rd1$delta_beta > 0, "Probe_ID"]])

        mean(fracs1, na.rm = TRUE)
    }, numeric(1)), decreasing = TRUE)
    sprintf("[%s](%1.1f) [%s](%1.1f)",
        names(fracs)[1], fracs[1], names(fracs)[2], fracs[2])
    
    ## results <- results[!(names(results) %in% ignore_branches)]
    ## cd <- meta[match(colnames(results), meta$betas),]
    ## se <- SummarizedExperiment(assays=list(results=results), colData=cd)
    ## metadata(se)$tissue_color <- metadata(reference)$tissue_color
    ## metadata(se)$branchID_color <- metadata(reference)$branchID_color
    ## se
}


