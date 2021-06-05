#' The openSesame pipeline
#'
#' This function is a simple wrapper of noob + nonlinear dye bias 
#' correction + pOOBAH masking.
#' 
#' If the input is an IDAT prefix or a \code{SigDF}, the output is
#' the beta value numerics. If the input is a minfi GenomicRatioSet
#' or RGChannelSet, the output is the sesamized GenomicRatioSet.
#' 
#' @param x SigDF(s), IDAT prefix(es), minfi GenomicRatioSet(s), 
#' or RGChannelSet(s)
#' @param platform optional platform string
#' @param manifest optional dynamic manifest
#' @param what either 'sigset' or 'beta'
#' @param ... parameters to getBetas
#' @param BPPARAM get parallel with MulticoreParam(n)
#' @return a numeric vector for processed beta values
#' @import BiocParallel
#' @examples
#' sdf <- sesameDataGet('HM450.1.TCGA.PAAD')$sdf
#' IDATprefixes <- searchIDATprefixes(
#'     system.file("extdata", "", package = "sesameData"))
#' betas <- openSesame(IDATprefixes)
#' @export
openSesame <- function(
    x, platform = '', manifest = NULL,
    BPPARAM=SerialParam(), ...) {

    ## expand if a directory
    if (length(x) == 1 && is(x, 'character') && dir.exists(x)) {
        x <- searchIDATprefixes(x)
    }
    
    if (is(x, "SigDF")) {
        getBetas(dyeBiasNL(noob(pOOBAH(qualityMask(x)))))
    } else if (is(x, "GenomicRatioSet")) {
        reopenSesame(x)
    } else if (is(x, "RGChannelSet")) {
        sesamize(x)
    } else if (is(x, 'character')) {
        if (length(x) == 1) {
            getBetas(dyeBiasNL(noob(pOOBAH(readIDATpair(
                x, platform = platform, manifest = manifest)
            ))))
        } else { # multiple IDAT prefixes / sigsets
            do.call(
                cbind, bplapply(x, openSesame,
                    platform = platform,
                    manifest = manifest, BPPARAM=BPPARAM, ...))
        }
    }
}

#' re-compute beta value for GenomicRatioSet
#' @param x GenomicRatioSet
#' @param naFrac  maximum NA fraction for a probe before it gets dropped (1)
#' @return a GenomicRatioSet
reopenSesame <- function(x, naFrac=0.2) { 
    pkgTest('minfi')
    stopifnot(is(x, "GenomicRatioSet"))
    if (!"Basename" %in% names(SummarizedExperiment::colData(x))) {
        stop("No column `Basename` in mcols(x)... cannot proceed.")
    } else { 
        message("Processing IDATs named in colData(x)$Basename...") 
    }
    
    SummarizedExperiment::assays(x)$Beta <- openSesame(
        x$Basename)[rownames(x),]
    keepRows <- names(which(rowSums(is.na(
        minfi::getBeta(x))) < round(ncol(x)*naFrac)))
    x[keepRows,]
}
