

#' Convert SNP from Infinium array to VCF file
#'
#' @param sset SigSet
#' @param vcf output VCF file path, if NULL output to console
#' @param refversion reference version, currently only support
#' @param annoS SNP variant annotation, download if not given
#' @param annoI Infinium-I variant annotation, download if not given
#' hg19 and hg38 in human
#'
#' @return VCF file. If vcf is NULL, a data.frame is output to
#' console. The data.frame does not contain VCF headers.
#' 
#' Note the vcf is not sorted. You can sort with
#' awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}'
#' 
#' @examples
#' sset <- sesameDataGet('EPIC.1.LNCaP')$sset
#'
#' annoS <- sesameDataPullVariantAnno_SNP('EPIC','hg19')
#' annoI <- sesameDataPullVariantAnno_InfiniumI('EPIC','hg19')
#' ## output to console
#' head(formatVCF(sset, annoS=annoS, annoI=annoI))
#' 
#' @export
formatVCF <- function(
    sset, vcf=NULL, refversion="hg19", annoS=NULL, annoI=NULL) {

    platform <- sset@platform

    if (is.null(annoS)) annoS <- sesameDataPullVariantAnno_SNP(platform)
    betas <- getBetas(sset, quality.mask=FALSE)[names(annoS)]
    vafs <- ifelse(annoS$U == 'REF', betas, 1-betas)
    vcflines_snp <- cbind(as.character(GenomicRanges::seqnames(annoS)),
        as.character(end(annoS)),
        names(annoS),
        annoS$REF, annoS$ALT,
        ".",
        ".",
        sprintf("VAF=%1.3f", vafs))


    if (is.null(annoI)) annoI <- sesameDataPullVariantAnno_InfiniumI(platform)
    af <- c(
        pmax(rowSums(oobR(sset)),1)/(
            pmax(rowSums(oobR(sset))+rowSums(IG(sset)),2)),
        pmax(rowSums(oobG(sset)),1)/(
            pmax(rowSums(oobG(sset))+rowSums(IR(sset)),2)))

    af <- af[names(annoI)]
    vafs <- ifelse(annoI$In.band == 'REF', af, 1-af)
    vcflines_typeI <- cbind(as.character(GenomicRanges::seqnames(annoI)),
        as.character(end(annoI)),
        names(annoI),
        annoI$REF, annoI$ALT,
        ".",
        ".",
        sprintf("PVF=%1.3f", vafs))

    header <- c(
        '##fileformat=VCFv4.0',
        sprintf('##fileDate=%s',format(Sys.time(),"%Y%m%d")),
        sprintf('##reference=%s', refversion),
        paste0(
            '##INFO=<ID=PVF,Number=1,Type=Float,',
            'Description="Pseudo Variant Frequency">'))
    
    out <- data.frame(rbind(
        vcflines_snp, vcflines_typeI))
    colnames(out) <- c(
        "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
    rownames(out) <- out$ID
    
    if(is.null(vcf)) {
        out
    } else {
        writeLines(header, vcf)
        write.table(
            out, file=vcf, append=TRUE, sep='\t',
            row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
}
