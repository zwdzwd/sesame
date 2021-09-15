
#' Check sample identity using SNP probes
#' 
#' @param betas numeric matrix (row: probes, column: samples)
#' @import wheatmap
#' @return grid object plotting SNP clustering
#' @examples 
#' betas <- sesameDataGet('HM450.10.TCGA.PAAD.normal')
#' SNPcheck(betas)
#' @export
SNPcheck <- function(betas) {
    probes <- names(sesameDataGet('HM450.probeInfo')$probe2chr.hg19)
    snp.probes <- probes[grep('rs',probes)]
    snp.probes <- intersect(snp.probes, rownames(betas))
    clus <- both.cluster(betas[snp.probes,])
    pkgTest('wheatmap')
    WHeatmap(
        clus$mat, 
        cmp=CMPar(dmin=0, dmax=1, stop.points = c('blue','white','red')),
        xticklabels=TRUE, yticklabels=TRUE, xticklabel.rotat = 45) + 
        WDendrogram(
            clus$column.clust, TopOf(height=0.2), facing='bottom')
}
