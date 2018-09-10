#' The openSesame pipeline
#'
#' This function is a simple wrapper of noob + nonlinear dye bias 
#' correction + pOOBAH masking.
#' 
#' If the input is an IDAT prefix or a \code{SigSet}, the output is
#' the beta value numerics. If the input is a minfi GenomicRatioSet
#' or RGChannelSet, the output is the sesamized GenomicRatioSet.
#' 
#' @param x SigSet(s), IDAT prefix(es), minfi GenomicRatioSet(s), 
#' or RGChannelSet(s)
#' @param ... parameters to getBetas
#' @return a numeric vector for processed beta values
#' @examples
#' sset <- sesameDataGet('HM450.1.TCGA.PAAD')$sset
#' IDATprefixes <- searchIDATprefixes(
#'     system.file("extdata", "", package = "sesameData"))
#' betas <- openSesame(IDATprefixes)
#' @export
openSesame <- function(x, ...) {
    if (length(x) == 1) {
        if (is(x, 'character')) {
            if (dir.exists(x)) {
                do.call(cbind, mclapply(searchIDATprefixes(x), openSesame))
            } else { # IDAT prefix
                x <- readIDATpair(x)
                stopifnot(is(x, 'SigSet'))
                getBetas(dyeBiasCorrTypeINorm(noob(x)), ...)
            }
        }
    } else if (is(x, "GenomicRatioSet")) {
        reopenSesame(x)
    } else if (is(x, "RGChannelSet")) {
        sesamize(x)
    } else { 
        do.call(cbind, mclapply(x, openSesame))
    }
}

#' re-compute beta value for GenomicRatioSet
#' @param x GenomicRatioSet
#' @param naFrac  maximum NA fraction for a probe before it gets dropped (1)
#' @return a GenomicRatioSet
reopenSesame <- function(x, naFrac=0.2) { 
    pkgTest('minfi')
    pkgTest('SummarizedExperiment')
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
