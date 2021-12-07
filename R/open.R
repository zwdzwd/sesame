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

    ## TODO add inferInfiniumIChannel
    if (is(x, "SigDF")) {
        getBetas(dyeBiasNL(noob(pOOBAH(qualityMask(x)))))
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
