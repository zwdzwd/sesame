#' The openSesame pipeline
#'
#' This function is a simple wrapper of noob + nonlinear dye bias 
#' correction + pOOBAH masking.
#' 
#' If the input is an IDAT prefix or a \code{SigDF}, the output is
#' the beta value numerics.
#' 
#' @param x SigDF(s), IDAT prefix(es)
#' @param platform optional platform string
#' @param manifest optional dynamic manifest
#' @param func either getBetas or getAFs
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
    x, platform = '', manifest = NULL, func = getBetas,
    BPPARAM=SerialParam(), ...) {

    ## expand if a directory
    if (length(x) == 1 && is(x, 'character') && dir.exists(x)) {
        x <- searchIDATprefixes(x)
    }

    ## TODO add inferInfiniumIChannel
    if (is(x, "SigDF")) {
        func(dyeBiasNL(noob(pOOBAH(qualityMask(x)))), ...)
    } else if (is(x, 'character')) {
        if (length(x) == 1) {
            func(dyeBiasNL(noob(pOOBAH(readIDATpair(
                x, platform = platform, manifest = manifest)
            ))))
        } else { # multiple IDAT prefixes / SigDFs
            do.call(
                cbind, bplapply(x, openSesame,
                    platform = platform, fun = func,
                    manifest = manifest, BPPARAM=BPPARAM, ...))
        }
    } else if (is(x, "list") && is(x[[1]], "SigDF")) {
        do.call(
            cbind, bplapply(x, openSesame,
                platform = platform, fun = func,
                manifest = manifest, BPPARAM=BPPARAM, ...))
    } else {
        stop("Unsupported input")
    }
}

