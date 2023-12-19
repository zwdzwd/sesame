## very simple genotyper
genotyper <- function(x, model_background=0.1, model_nbeads=40) {

    GL <- vapply(
        c(model_background, 0.5, 1-model_background),
        function(af) {
            dbinom(
                round(x*model_nbeads),
                size=model_nbeads, prob=af)}, numeric(1))
        
    ind <- which.max(GL)
    GT <- c('0/0','0/1','1/1')[ind]
    GS <- floor(-log10(1-GL[ind] / sum(GL))*10) # assuming equal prior
    list(GT=GT, GS=GS)
}

vcf_header <- function(genome) {
    c('##fileformat=VCFv4.0',
        sprintf('##fileDate=%s',format(Sys.time(),"%Y%m%d")),
        sprintf('##reference=%s', genome),
        paste0('##INFO=<ID=PVF,Number=1,Type=Float,',
            'Description="Pseudo Variant Frequency">'),
        paste0('##INFO=<ID=GT,Number=1,Type=String,',
            'Description="Genotype">'),
        paste0('##INFO=<ID=GS,Number=1,Type=Integer,',
            'Description="Genotyping score from 7 to 85">'),
        paste0('##INFO=<ID=Probe_ID,Number=1,Type=String,',
            'Description="Infinium Probe ID">'),
        paste0('##INFO=<ID=rs_ID,Number=1,Type=String,',
            'Description="Overlapping rs ID from dbSNP">'))
}

#' Convert SNP from Infinium array to VCF file
#'
#' @param sdf SigDF
#' @param anno SNP variant annotation, available at
#' https://github.com/zhou-lab/InfiniumAnnotationV1/tree/main/Anno/EPIC
#' EPIC.hg38.snp.tsv.gz
#' @param vcf output VCF file path, if NULL output to console
#' @param genome genome
#' @param verbose print more messages
#' @return VCF file. If vcf is NULL, a data.frame is output to
#' console. The data.frame does not contain VCF headers.
#' Note the output vcf is not sorted.
#' 
#' @importFrom utils write.table
#' @examples
#' sesameDataCacheAll() # if not done yet
#' sdf <- sesameDataGet('EPIC.1.SigDF')
#'
#' \dontrun{
#' ## download anno from
#' ## http://zwdzwd.github.io/InfiniumAnnotation
#' ## output to console
#' anno = read_tsv(sesameAnno_download("EPICv2.hg38.snp.tsv.gz"))
#' head(formatVCF(sdf, anno))
#' }
#' 
#' @export
formatVCF <- function(
    sdf, anno, vcf=NULL, genome="hg38", verbose = FALSE) {
    
    platform <- sdfPlatform(sdf, verbose = verbose)
    betas <- getBetas(sdf)[anno$Probe_ID]
    af <- getAFTypeIbySumAlleles(sdf, known.ccs.only=FALSE)
    vafs <- ifelse(anno$U == "ALT", 1-betas, betas)
    vafs <- ifelse(anno$U == "REF_InfI", af[anno$Probe_ID], vafs)
    
    gts <- lapply(vafs, genotyper)
    GT <- vapply(gts, function(g) g$GT, character(1))
    GS <- vapply(gts, function(g) g$GS, numeric(1))
    anno$REF[anno$REF == "ACT"] <- "H"
    anno$REF[anno$REF == "AGT"] <- "D"
    anno$ALT[anno$ALT == "ACT"] <- "H"
    anno$ALT[anno$ALT == "AGT"] <- "D"
    vcflines <- cbind(anno$chrm, anno$end,
        ".", anno$REF, anno$ALT, GS, ifelse(GS>20,'PASS','FAIL'),
        paste0(sprintf(
            "PVF=%1.3f;GT=%s;GS=%d;Probe_ID=%s",
            vafs, GT, GS, anno$Probe_ID),
        ifelse(is.na(anno$rs), "", paste0(";rs_ID=", anno$rs))))

    header <- vcf_header(genome)    
    out <- data.frame(vcflines)
    colnames(out) <- c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
    out <- out[order(out[['#CHROM']], as.numeric(out[['POS']])),]
    
    if(is.null(vcf)) { return(out);
    } else {
        writeLines(header, vcf)
        write.table(out, file=vcf, append=TRUE, sep='\t',
            row.names = FALSE, col.names = FALSE, quote = FALSE) }
}
