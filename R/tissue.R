
#' Compare mouse array data with mouse tissue references
#'
#' @param betas matrix of betas for the target sample
#' @return grid object that contrast the target sample with
#' pre-built mouse tissue reference
#' @import SummarizedExperiment
#' @export
#' @examples
#' b = sesameDataGet("MM285.10.tissue")$betas[,1:2]
#' compareMouseTissueReference(b)
compareMouseTissueReference = function(betas) {
    se = sesameDataGet("MM285.tissueSignature")
    cd = as_tibble(colData(se))
    rd = as_tibble(rowData(se))
    md = metadata(se)

    if (is.null(dim(betas))) {
        betas = cbind(betas)
    }
    
    WHeatmap(assay(se), cmp=CMPar(stop.points=c('blue','yellow')),
        name="b1") + 
        WHeatmap(betas[rownames(se),], RightOf("b1"),
            cmp=CMPar(stop.points=c('blue','yellow')), name="b2") +
        WColorBarH(cd$Tissue_Corrected, TopOf("b1",height=0.03),
            cmp=CMPar(label2color=md$tissue_color), name="ti") +
        WColorBarV(rd$probebranches, RightOf("b2", width=0.1),
            cmp=CMPar(label2color=md$branch_color), name="bh") +
        WLegendV("ti", TopRightOf("bh", just=c('left','top'), h.pad=0.02),
            height=0.02) + 
        WLegendV('bh', Beneath(pad=0.06)) + 
        WCustomize(mar.bottom=0.15, mar.right=0.06)
}
