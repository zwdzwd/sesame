
#' Apply a chain of sesame preprocessing functions in an arbitrary order
#' 
#' @param sdf SigDF
#' @param spell code that indicates preprocessing functions and their
#' execution order (functions on the left is executed first).
#' @return SigDF
#' @examples
#' sdf <- sesameDataGet("MM285.1.SigDF")
#' sdf1 <- prepSesame(sdf, "QCDPB")
#' @export
prepSesame <- function(sdf, spell = "QCDPB") {
    sp2fun <- list(
        Q = qualityMask,
        C = inferInfiniumIChannel,
        D = dyeBiasNL,
        P = pOOBAH,
        B = noob,
        S = inferSpecies,
        T = inferStrain,
        M = matchDesign)
    
    spells <- str_split(spell,"")[[1]]
    stopifnot(all(spells %in% names(sp2fun)))
    x <- sdf
    for(sp1 in spells) {
        x <- sp2fun[[sp1]](x) }
    x
}

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
#' @param preprocess preprocessing spell, see preprocessSesame()
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
    x, platform = '', preprocess = "QCDPB", manifest = NULL,
    func = getBetas, BPPARAM=SerialParam(), ...) {

    ## expand if a directory
    if (length(x) == 1 && is(x, 'character') && dir.exists(x)) {
        x <- searchIDATprefixes(x)
    }

    if (is(x, "SigDF")) {
        func(prepSesame(x), ...)
    } else if (is(x, 'character')) {
        if (length(x) == 1) {
            func(prepSesame(readIDATpair(
                x, platform = platform, manifest = manifest)), ...)
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

