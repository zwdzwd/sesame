reference_plot_se = function(betas, se, color=c("blueYellow","fullJet")) {

    pkgTest("wheatmap")
    color = match.arg(color)
    if (color == "blueYellow") stop.points = c("blue","yellow")
    else stop.points = NULL

    cd = as_tibble(colData(se))
    rd = as_tibble(rowData(se))
    md = metadata(se)
    if (!is.null(betas) && is.null(dim(betas))) { # in case a vector
        betas = cbind(betas)
    }

    g = WHeatmap(assay(se), cmp=CMPar(stop.points=stop.points,
        dmin=0, dmax=1), name="b1") # reference
    if (!is.null(betas)) {          # query samples
        g = g + WHeatmap(betas[rownames(se),], RightOf("b1"),
            cmp=CMPar(stop.points=stop.points, dmin=0, dmax=1),
            name="b2", xticklabels=TRUE, xticklabels.n=ncol(betas))
        right = "b2"
    } else { # in case target is not given, plot just the reference
        right = "b1"
    }
    ## branch color bar (vertical)
    g = g + WColorBarV(rd$probebranches, RightOf(right, width=0.03),
        cmp=CMPar(label2color=md$branch_color), name="bh")
    ## tissue color bar (horizontal)
    g = g + WColorBarH(cd$Tissue_Corrected, TopOf("b1",height=0.03),
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
#' b = sesameDataGet("MM285.10.tissue")$betas[,1:2]
#' compareMouseTissueReference(b)
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
compareMouseTissueReference = function(betas=NULL, color="blueYellow") {
    se = sesameDataGet("MM285.tissueSignature")
    reference_plot_se(betas, se, color=color)
}

#' Compare beta value against mouse blood reference
#'
#' @param betas matrix of betas for the target sample
#' @param color either blueYellow or fullJet
#' @return grid object that co-plots with a pre-built mouse blood reference
#' @export
#' @examples
#' b = sesameDataGet("MM285.10.tissue")$betas[,10]
#' compareMouseBloodReference(b)
compareMouseBloodReference = function(betas=NULL, color="blueYellow") {
    se = sesameDataGet("MM285.bloodSignature")
    reference_plot_se(betas, se, color=color)
}

