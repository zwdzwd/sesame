
#' Turn beta values into a UCSC browser track
#'
#' @param betas a named numeric vector
#' @param output output file name
#' @param platform HM450, EPIC etc.
#' @param refversion hg38, hg19 etc.
#' @return when output is null, return a data.frame, otherwise NULL
#' @examples
#'
#' betas.tissue <- sesameDataGet('HM450.1.TCGA.PAAD')$betas
#' ## add output to create an actual file
#' df <- createUCSCtrack(betas.tissue)
#' 
#' ## to convert to bigBed
#' ## sort -k1,1 -k2,2n output.bed >output_sorted.bed
#' ## bedToBigBed output_sorted.bed hg38.chrom output.bb
#' @export
createUCSCtrack <- function(
    betas, output=NULL, platform='HM450', refversion='hg38') {
    
    probeInfo <- sesameDataGet(paste0(
        platform, '.', refversion, '.manifest'))

    probeInfo <- probeInfo[GenomicRanges::seqnames(probeInfo) != "*"]

    betas <- betas[names(probeInfo)]
    df <- data.frame(
        chrm = GenomicRanges::seqnames(probeInfo),
        beg = GenomicRanges::start(probeInfo)-1,
        end = GenomicRanges::end(probeInfo),
        name = names(probeInfo),
        score = ifelse(is.na(betas), 0, as.integer(betas*1000)),
        strand = GenomicRanges::strand(probeInfo),
        thickStart = GenomicRanges::start(probeInfo)-1,
        thickEnd = GenomicRanges::end(probeInfo),
        itemRgb = ifelse(
            is.na(betas), '0,0,0',
            ifelse(
                betas < 0.3, '0,0,255', # blue
                ifelse(
                    betas > 0.7, '255,0,0', # red
                    '50,150,0'))) # green
    )

    if (is.null(output))
        df
    else
        write.table(
            df, file=output, col.names=FALSE,
            row.names=FALSE, quote=FALSE, sep='\t')
}
