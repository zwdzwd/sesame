#' The openSesame pipeline
#'
#' This function is a simple wrapper of noob + nonlinear dye bias 
#' correction + pOOBAH masking.
#' 
#' If the input is an IDAT prefix or a \code{SigDF}, the output is
#' the beta value numerics.
#'
#' Notes on the order of operation:
#' 1. qualityMask and inferSpecies should go before noob and pOOBAH,
#' otherwise the background is too high because of Multi,
#' uk and other probes
#' 2. dyeBias correction needs to happen early
#' 3. channel inference before dyebias
#' 4. noob should happen last, pOOBAH before noob because noob modifies oob
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
        func(noob(pOOBAH(dyeBiasNL(inferInfiniumIChannel(qualityMask(x))))),
            ...)
    } else if (is(x, 'character')) {
        if (length(x) == 1) {
            func(noob(pOOBAH(dyeBiasNL(inferInfiniumIChannel(qualityMask(
                readIDATpair(x, platform = platform, manifest = manifest)
            ))))), ...)
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

