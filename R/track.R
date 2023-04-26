
#' Turn beta values into a UCSC browser track
#'
#' @param betas a named numeric vector
#' @param output output file name
#' @param platform HM450, EPIC etc.
#' @param genome hg38, mm10, ..., will infer if not given.
#' For additional mapping, download the GRanges object from
#' http://zwdzwd.github.io/InfiniumAnnotation
#' and provide the following argument
#' ..., genome = sesameAnno_buildManifestGRanges("downloaded_file"),...
#' to this function.
#' @return when output is null, return a data.frame, otherwise NULL
#' @importFrom utils write.table
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
    betas, output=NULL, platform='HM450', genome='hg38') {
    
    probeInfo <- sesameData_getManifestGRanges(platform, genome)

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
